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
#include <eigendecompose.h>

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

void Method::reset_mass_nodes(Grid &grid, int in)
{
  if (should_compute_mass_nodes())
    grid.mass.at(in) = 0;
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
    double normalized_wf = wf*solid.mass.at(ip)/node_mass;
    const Vector3d &dx = solid.grid->x0.at(in) - solid.x.at(ip);
    
    Vector3d vtemp = solid.v.at(ip);

    if (apic() || update->method->ge)
    {
      if (is_TL)
        vtemp += solid.Fdot.at(ip)*(solid.grid->x0.at(in) - solid.x0.at(ip));
      else
        vtemp += solid.L.at(ip)*(solid.grid->x0.at(in) - solid.x.at(ip));
    }

    solid.grid->v.at(in) += normalized_wf*vtemp;

    if (solid.grid->rigid.at(in))
      solid.grid->mb.at(in) += normalized_wf*solid.v_update.at(ip);

    if (temp)
      solid.grid->T.at(in) += normalized_wf*solid.T.at(ip);
  }
}

void Method::compute_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd)
{
  compute_internal_force_nodes(solid, in, ip, wf, wfd);

  if (!solid.grid->rigid.at(in))
    solid.grid->mb.at(in) += wf*solid.mbp.at(ip);

  if (temp)
  {
    if (solid.domain->axisymmetric)
    {
      error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
    }

    if (solid.grid->mass.at(in))
      solid.grid->Qext.at(in) += wf*solid.gamma.at(ip);

    solid.grid->Qint.at(in) += wfd.dot(solid.q.at(ip));
  }
}

void Method::update_grid_velocities(Grid &grid, int in)
{
  double T_update;

  Vector3d &v_update = grid.v_update.at(in) = grid.v.at(in);
  if (temp)
    T_update = grid.T_update.at(in) = grid.T.at(in);

  if (double mass = grid.mass.at(in))
  {
    if (!grid.rigid.at(in))
      v_update += update->dt*(grid.f.at(in) + grid.mb.at(in))/mass;

    if (temp)
      T_update += update->dt*(grid.Qint.at(in) + grid.Qext.at(in))/mass;
  }
}

void Method::compute_velocity_acceleration(Solid &solid, int in, int ip, double wf)
{
  solid.v_update.at(ip) += wf*solid.grid->v_update.at(in);

  if (solid.mat->rigid)
    return;

  const Vector3d &delta_a = wf*(solid.grid->v_update.at(in) - solid.grid->v.at(in))/update->dt;
  solid.a.at(ip) += delta_a;
  solid.f.at(ip) += delta_a*solid.mass.at(ip);

  if (temp)
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
    const Vector3d &dx = is_TL? solid.grid->x0.at(in) - solid.x0.at(ip):
                                solid.grid->x.at(in) - solid.grid->x0.at(in);

    Matrix3d gradient;
    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradient(j, k) += vn.at(in)[j]*dx[k]*wf;

    gradients.at(ip) += gradient*solid.Di;
  }
  else
    for (int j = 0; j < domain->dimension; j++)
      for (int k = 0; k < domain->dimension; k++)
        gradients.at(ip)(j, k) += vn.at(in)[j]*wfd[k];

  if (domain->dimension == 2 && domain->axisymmetric)
    gradients.at(ip)(2, 2) += vn.at(in)[0]*wf/solid.x0.at(ip)[0];

  if (temp)
    solid.q.at(ip) -= wfd*(doublemapping? solid.grid->T: solid.grid->T_update).at(in)
                      *(is_TL? solid.vol0: solid.vol).at(ip)*solid.mat->invcp*solid.mat->kappa;
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

// EB: remove this
#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333

void Method::update_stress(bool doublemapping, Solid &solid, int ip)
{
  if (solid.mat->rigid)
    return;

  if (solid.mat->type == material->constitutive_model::LINEAR)
  {
    const Matrix3d &strain_increment = update->dt*solid.D.at(ip);
    solid.strain_el.at(ip) += strain_increment;
    solid.sigma.at(ip) += 2*solid.mat->G*strain_increment +
      solid.mat->lambda*strain_increment.trace()*Matrix3d::identity();

    if (is_TL)
      solid.vol0PK1.at(ip) = solid.vol0.at(ip)*solid.J.at(ip)*solid.R.at(ip)*solid.sigma.at(ip)*
                            solid.R.at(ip).transpose()*solid.Finv.at(ip).transpose();
    solid.gamma.at(ip) = 0;
  }
  else if (solid.mat->type == material->constitutive_model::NEO_HOOKEAN)
  {
    // Neo-Hookean material:
    const Matrix3d &FinvT = solid.Finv.at(ip).transpose();
    const Matrix3d &PK1 = solid.mat->G*(solid.F.at(ip) - FinvT) + solid.mat->lambda*log(solid.J.at(ip))*FinvT;
    solid.vol0PK1.at(ip) = solid.vol0.at(ip)*PK1;
    solid.sigma.at(ip) = 1/solid.J.at(ip)*solid.F.at(ip)*PK1.transpose();

    solid.strain_el.at(ip) = 0.5*(solid.F.at(ip).transpose()*solid.F.at(ip) - Matrix3d::identity());
    solid.gamma.at(ip) = 0;
  }
  else
  {
    double pH = 0;
    double plastic_strain_increment = 0;
    Matrix3d sigma_dev;

    double T = solid.mat->cp? solid.T.at(ip): 0;

    solid.mat->eos->compute_pressure(pH, solid.ienergy.at(ip), solid.J.at(ip), solid.rho.at(ip),
                                      solid.damage.at(ip), solid.D.at(ip), solid.grid->cellsize, T);

    if (solid.mat->cp)
      pH += solid.mat->alpha*(solid.T.at(ip) - solid.T0);

    sigma_dev = solid.mat->strength->update_deviatoric_stress(
      solid.sigma.at(ip), solid.D.at(ip), plastic_strain_increment,
      solid.eff_plastic_strain.at(ip), solid.eff_plastic_strain_rate.at(ip), solid.damage.at(ip),
      T);

    solid.eff_plastic_strain.at(ip) += plastic_strain_increment;

    // compute a characteristic time over which to average the plastic strain
    solid.eff_plastic_strain_rate.at(ip) += (plastic_strain_increment - solid.eff_plastic_strain_rate.at(ip)*update->dt)/
                                            1000/solid.grid->cellsize*solid.mat->signal_velocity;
    solid.eff_plastic_strain_rate.at(ip) = MAX(0.0, solid.eff_plastic_strain_rate.at(ip));

    if (solid.mat->damage)
        solid.mat->damage->compute_damage(solid.damage_init.at(ip), solid.damage.at(ip), pH,
                                          sigma_dev, solid.eff_plastic_strain_rate.at(ip),
                                          plastic_strain_increment, T);

    if (solid.mat->temp)
    {
      solid.mat->temp->compute_heat_source(solid.T.at(ip), solid.gamma.at(ip), SQRT_3_OVER_2*sigma_dev.norm(),
                                           solid.eff_plastic_strain_rate.at(ip));
      if (is_TL)
        solid.gamma.at(ip) *= solid.vol0.at(ip)*solid.mat->invcp;
      else
        solid.gamma.at(ip) *= solid.vol.at(ip)*solid.mat->invcp;
    }
    else
	  solid.gamma.at(ip) = 0;

    solid.sigma.at(ip) = -pH*(1 - (pH < 0? solid.damage.at(ip): 0))*Matrix3d::identity() + sigma_dev;

      solid.strain_el.at(ip) =
        (update->dt*solid.D.at(ip).trace() + solid.strain_el.at(ip).trace())/3*Matrix3d::identity() +
        sigma_dev/solid.mat->G/(1 - (solid.damage.at(ip) > 1e-10? solid.damage.at(ip): 0));

    if (is_TL)
      solid.vol0PK1.at(ip) = solid.vol0.at(ip)*solid.J.at(ip)*
        solid.R.at(ip)*solid.sigma.at(ip)*solid.R.at(ip).transpose()*
        solid.Finv.at(ip).transpose();
  }

  if (solid.damage.at(ip) >= 1.0)
    return;

  double p_wave_speed = sqrt((solid.mat->K + FOUR_THIRD*solid.mat->G)/solid.rho.at(ip)) +
                        MAX(MAX(fabs(solid.v.at(ip)(0)), fabs(solid.v.at(ip)(1))), fabs(solid.v.at(ip)(2)));

  if (std::isnan(p_wave_speed))
  {
    cout << "Error: max_p_wave_speed is nan with ip=" << ip
      << ", ptag.at(ip)=" << solid.ptag.at(ip) << ", rho0.at(ip)=" << solid.rho0.at(ip)<< ", rho.at(ip)=" << solid.rho.at(ip)
      << ", K=" << solid.mat->K << ", G=" << solid.mat->G << ", J.at(ip)=" << solid.J.at(ip)
      << endl;
    error->one(FLERR, "");
  }
  else if (p_wave_speed < 0.0)
  {
    cout << "Error: p_wave_speed= " << p_wave_speed
      << " with ip=" << ip << ", rho.at(ip)=" << solid.rho.at(ip) << ", K=" << solid.mat->K
      << ", G=" << solid.mat->G << endl;
    error->one(FLERR, "");
  }
  
  double h_ratio = 1;

  if (is_TL)
  {
    Matrix3d eigenvalues = solid.F.at(ip);
    eigendecompose(eigenvalues);
    // EB: revisit
    if (/*esF.info()== Success*/false)
    {
      h_ratio = MIN(h_ratio, fabs(eigenvalues(0, 0)));
      h_ratio = MIN(h_ratio, fabs(eigenvalues(1, 1)));
      h_ratio = MIN(h_ratio, fabs(eigenvalues(2, 2)));
    }

    if (h_ratio == 0)
    {
      cout << "min_h_ratio == 0 with ip=" << ip
        << "F=\n" <<  solid.F.at(ip) << endl
        << "eigenvalues of F:" << eigenvalues(0, 0) << "\t" << eigenvalues(1, 1) << "\t" << eigenvalues(2, 2) << endl;
      //cout << "esF.info()=" << esF.info() << endl;
      error->one(FLERR, "");
    }
  }
  else
  {
    h_ratio = 1;
  }

  solid.dtCFL.at(ip) = solid.grid->cellsize*h_ratio/p_wave_speed;

  if (std::isnan(solid.dtCFL.at(ip)))
  {
    cout << "Error: dtCFL = " << solid.dtCFL.at(ip) << "\n";
    cout << "p_wave_speed = " << p_wave_speed
      << ", grid->cellsize=" << solid.grid->cellsize << endl;
    error->one(FLERR, "");
  }
}

void Method::adjust_dt()
{
  if (update->dt_constant) return; // dt is set as a constant, do not update

  double dtCFL = 1.0e22;
  double dtCFL_reduced = 1.0e22;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
    for (int ip = 0; ip < domain->solids.at(isolid)->np_local; ip++)
    {
      dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL.at(ip));
      if (!dtCFL)
      {
        cout << "Error: dtCFL == 0\n";
        cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
        error->one(FLERR, "");
      }
      else if (std::isnan(dtCFL))
      {
        cout << "Error: dtCFL = " << dtCFL << "\n";
        cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL.at(ip) << "\n";
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
    for (int ip = 0; ip < domain->solids[isolid]->np_local; ip++)
      domain->solids[isolid]->mbp[ip] = Vector3d();
  }
}
