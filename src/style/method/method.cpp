/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <method.h>

#include <solid.h>
#include <update.h>
#include <domain.h>
#include <error.h>
#include <universe.h>
#include <input.h>
#include <inverse.h>

using namespace std;

Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
  ge = false;
  temp = false;
}

bool Method::apic()
{
  switch (update->sub_method_type)
  {
  case Update::SubMethodType::APIC:
  case Update::SubMethodType::AFLIP:
  case Update::SubMethodType::ASFLIP:
  case Update::SubMethodType::MLS:
    return true;
  }

  return false;
}

void Method::compute_mass_nodes(Solid &solid, int in, int ip, double wf)
{
  if (should_compute_mass_nodes() && !solid.grid->rigid.at(in) || solid.mat->rigid)
    solid.grid->mass.at(in) += wf*solid.mass.at(ip);
}

void Method::compute_velocity_nodes(Solid &solid, int in, int ip, double wf)
{
  if (solid.grid->rigid.at(in) && !solid.mat->rigid)
    return;

  if (double node_mass = solid.grid->mass.at(in))
  {
    double wf_mass = wf*solid.mass.at(ip);
    const Vector3d &dx = solid.grid->x0.at(in) - solid.x.at(ip);
    
    Vector3d vtemp = solid.v.at(ip);

    if (apic())
    {
      if (is_TL)
        vtemp += solid.Fdot.at(ip)*(solid.grid->x0.at(in) - solid.x0.at(ip));
      else
        vtemp += solid.L.at(ip)*dx;

      solid.grid->v.at(in) += vtemp/node_mass;
    }
    else
    {
      if (update->method->ge)
        vtemp += solid.L.at(ip)*dx;

      solid.grid->v.at(in) += wf_mass*vtemp/node_mass;
    }

    if (solid.grid->rigid.at(in))
      solid.grid->mb.at(in) += wf_mass*solid.v_update.at(ip)/node_mass;
  }
}

void Method::compute_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  compute_internal_force_nodes(solid, in, ip, wf, wfd);
  //if (TL)
  //{
  //  grid->f.at(in) -= vol0PK1.at(ip)*wfd;

  //  if (domain->axisymmetric)
  //    grid->f.at(in)[0] -= vol0PK1.at(ip)(2, 2)*wf/x0.at(ip)[0];
  //}
  //else
  //{
  //  if (MLS)
  //    grid->f.at(in) -= vol.at(ip)*wf*sigma.at(ip)*Di*
  //                      (grid->x0.at(in) - (is_TL? x0: x).at(ip));
  //  else
  //    grid->f.at(in) -= vol.at(ip)*sigma.at(ip)*wfd;

  //  if (domain->axisymmetric)
  //    grid->f.at(in)[0] -= vol.at(ip)*sigma.at(ip)(2, 2)*wf/x.at(ip)[0];
  //}

  if (!solid.grid->rigid.at(in))
    solid.grid->mb.at(in) += wf*solid.mbp.at(ip);
}

void Method::compute_temperature_nodes(Solid &solid, int in, int ip, double wf)
{
  if (double node_mass = solid.grid->mass.at(in))
    solid.grid->T.at(in) += wf*solid.mass.at(ip)*solid.T.at(ip)/node_mass;
}

void Method::compute_temperature_driving_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  if (solid.domain->axisymmetric)
  {
    error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
  }

  if (solid.grid->mass.at(in))
    solid.grid->Qext.at(in) += wf*solid.gamma.at(ip);

  solid.grid->Qint.at(in) += wfd.dot(solid.q.at(ip));
}

void Method::update_grid_velocities(Grid &grid, int in)
{
  Vector3d &v_update = grid.v_update.at(in) = grid.v.at(in);

  if (!grid.rigid.at(in))
    if (double mass = grid.mass.at(in))
      v_update += update->dt*(grid.f.at(in) + grid.mb.at(in))/mass;
}

void Method::update_grid_temperature(Grid &grid, int in)
{
  double T_update = grid.T_update.at(in) = grid.T.at(in);

  if (double mass = grid.mass.at(in))
     T_update += update->dt*(grid.Qint.at(in) + grid.Qext.at(in))/mass;
}

void Method::compute_velocity_acceleration(Solid &solid, int in, int ip, double wf)
{
  solid.v_update.at(ip) += wf*solid.grid->v_update.at(in);

  if (solid.mat->rigid)
    return;

  const Vector3d &delta_a = wf*(solid.grid->v_update.at(in) - solid.grid->v.at(in))/update->dt;
  solid.a.at(ip) += delta_a;
  solid.f.at(ip) += delta_a*solid.mass.at(ip);
}

void Method::compute_particle_temperature(Solid &solid, int in, int ip, double wf)
{
  solid.T.at(ip) += wf*solid.grid->T_update.at(in);
}

void Method::update_position(Solid &solid, int ip)
{
  check_particle_in_domain(solid.x.at(ip) += update->dt*solid.v_update.at(ip), ip);
}

void Method::advance_particles(Solid &solid, int ip)
{
  Vector3d &v = solid.v.at(ip);

  if (update->sub_method_type == Update::SubMethodType::ASFLIP)
    solid.x.at(ip) += update->dt*v;

  v = (1 - update->PIC_FLIP)*solid.v_update.at(ip) + update->PIC_FLIP*(v + update->dt*solid.a.at(ip));
}

      //solid.reset_rate_deformation_gradient(false);
void Method::compute_rate_deformation_gradient(bool doublemapping, Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  if (solid.mat->rigid)
    return;
        
  vector<Matrix3d> &gradients = get_gradients(solid);
  const vector<Vector3d> &vn = doublemapping? solid.grid->v: solid.grid->v_update;

  if (update->sub_method_type == Update::SubMethodType::APIC)
  {
    const Vector3d &dx = solid.grid->x.at(in) - solid.grid->x0.at(in);

    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradients.at(ip)(j, k) += vn.at(in)[j]*dx[k]*wf;
  }
  else
    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradients.at(ip)(j, k) += vn.at(in)[j]*wfd[k];

  if (domain->dimension == 2 && domain->axisymmetric)
    gradients.at(ip)(2, 2) += vn.at(in)[0]*wf/solid.x0.at(ip)[0];
}

void Method::update_deformation_gradient_matrix(Solid &solid, int ip)
{
  // FOR CPDI:
  //if (update->method->style == 1)
  //  solid.F.at(ip) += update->dt*solid.Fdot.at(ip);
  //else

  solid.F.at(ip) = (Matrix3d::identity() + update->dt*solid.L.at(ip))*solid.F.at(ip);
}

void Method::update_deformation_gradient_determinant(Solid &solid, int ip)
{
  // FOR CPDI:
  //if (update->method->style == 1)
  //{
  //  solid.vol.at(ip) = 0.5*(solid.xpc[solid.nc*ip + 0][0]*solid.xpc[solid.nc*ip + 1][1] -
  //                          solid.xpc[solid.nc*ip + 1][0]*solid.xpc[solid.nc*ip + 0][1] +
  //                          solid.xpc[solid.nc*ip + 1][0]*solid.xpc[solid.nc*ip + 2][1] -
  //                          solid.xpc[solid.nc*ip + 2][0]*solid.xpc[solid.nc*ip + 1][1] +
  //                          solid.xpc[solid.nc*ip + 2][0]*solid.xpc[solid.nc*ip + 3][1] -
  //                          solid.xpc[solid.nc*ip + 3][0]*solid.xpc[solid.nc*ip + 2][1] +
  //                          solid.xpc[solid.nc*ip + 3][0]*solid.xpc[solid.nc*ip + 0][1] -
  //                          solid.xpc[solid.nc*ip + 0][0]*solid.xpc[solid.nc*ip + 3][1]);
  //  solid.J.at(ip) = solid.vol.at(ip)/solid.vol0.at(ip);
  //}
  //else
  //{

  solid.J.at(ip) = determinant(solid.F.at(ip));
  solid.vol.at(ip) = solid.J.at(ip)*solid.vol0.at(ip);
}

void Method::update_deformation_gradient(Solid &solid, int ip)
{
  if (solid.mat->rigid)
    return;

  update_deformation_gradient_matrix(solid, ip);

  solid.Finv.at(ip) = inverse(solid.F.at(ip));

  update_deformation_gradient_determinant(solid, ip);

  if (solid.J.at(ip) <= 0.0 && solid.damage.at(ip) < 1.0)
  {
    cout << "Error: J[" << solid.ptag.at(ip) << "]<=0.0 == " << solid.J.at(ip) << endl;
    cout << "F[" << solid.ptag.at(ip) << "]:" << endl << solid.F.at(ip) << endl;
    cout << "Fdot[" << solid.ptag.at(ip) << "]:" << endl << solid.Fdot.at(ip) << endl;
    cout << "damage[" << solid.ptag.at(ip) << "]:" << endl << solid.damage.at(ip) << endl;
    error->one(FLERR, "");
  }

  solid.rho.at(ip) = solid.rho0.at(ip)/solid.J.at(ip);

  if (solid.mat->type != material->constitutive_model::NEO_HOOKEAN)
    update_velocity_gradient_matrix(solid, ip);
}

void Method::update_stress(bool doublemapping)
{
  for (Solid *solid: domain->solids)
  {
    solid->update_stress();

    if (temp)
    {
      solid->reset_heat_flux();

      for (int i = 0; i < solid->neigh_n.size(); i++)
      {
        int in = solid->neigh_n.at(i);
        int ip = solid->neigh_p.at(i);
        const Vector3d &wfd = solid->wfd.at(i);

        solid->compute_heat_flux(in, ip, wfd, doublemapping);
      }
    }
  }
}

void Method::adjust_dt()
{
  if (update->dt_constant) return; // dt is set as a constant, do not update

  double dtCFL = 1.0e22;
  double dtCFL_reduced = 1.0e22;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0)
    {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      error->one(FLERR, "");
    } else if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      error->one(FLERR, "");
    }
  }

  MPI_Allreduce(&dtCFL, &dtCFL_reduced, 1, MPI_DOUBLE, MPI_MIN, universe->uworld);

  update->dt = dtCFL_reduced * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void Method::reset()
{
  int np_local;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np_local = domain->solids[isolid]->np_local;
    for (int ip = 0; ip < np_local; ip++) domain->solids[isolid]->mbp[ip] = Vector3d();
  }
}
