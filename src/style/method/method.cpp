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
#include <basis_functions.h>

using namespace std;

Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
  ge = false;
  temp = false;
}

void Method::compute_grid_weight_functions_and_gradients(Solid &solid, int ip)
{
  if (update->ntimestep == 0 && solid.mat->rigid)
    rigid_solids = 1;

  bigint nnodes = solid.grid->nnodes;
  bigint nnodes_local = solid.grid->nnodes_local;
  bigint nnodes_ghost = solid.grid->nnodes_ghost;

  Vector3d r;
  double s[3], sd[3];
  const Vector3d &xp = solid.x[ip];
  Kokkos::View<Vector3d*> &x0 = solid.grid->x0;
  double inv_cellsize = 1/solid.grid->cellsize;
  double wf;
  Vector3d wfd;

  Kokkos::View<Vector3i*> &ntype = solid.grid->ntype;

  vector<tagint> &map_ntag = solid.grid->map_ntag;

  vector<int> n_neigh;
    
  int ny = solid.grid->ny_global;
  int nz = solid.grid->nz_global;

  if (nnodes_local + nnodes_ghost)
  {
    // Calculate what nodes particle ip will interact with:
    n_neigh.clear();

    int half_support = update->shape_function == Update::ShapeFunctions::LINEAR? 1: 2;

    int i0 = (xp[0] - domain->boxlo[0])*inv_cellsize + 1 - half_support;
    for (int i = i0; i < i0 + 2*half_support; i++)
    {
      if (ny > 1)
      {
        int j0 = (xp[1] - domain->boxlo[1])*inv_cellsize + 1 - half_support;
        for (int j = j0; j < j0 + 2*half_support; j++)
        {
          if (nz > 1)
          {
            int k0 = (xp[2] - domain->boxlo[2])*inv_cellsize + 1 - half_support;
            for (int k = k0; k < k0 + 2*half_support; k++)
            {
              int tag = nz*ny*i + nz*j + k;

              if (tag < nnodes)
              {
                tagint inn = map_ntag[tag];

                if (inn != -1)
                  n_neigh.push_back(inn);
              }
            }
          }
          else
          {
            int tag = ny*i + j;

            if (tag < nnodes)
            {
              tagint in = map_ntag[tag];

              if (in != -1)
                n_neigh.push_back(in);
            }
          }
        }
      }
      else
      {
        if (i < nnodes_local + nnodes_ghost)
          n_neigh.push_back(i);
      }
    }

    for (int i = 0; i < n_neigh.size(); i++)
    {
      int in = n_neigh[i];

      // Calculate the distance between each pair of particle/node:
      r = (xp - x0[in])*inv_cellsize;

      s[0] = basis_function(r[0], ntype[in][0]);
      wf = s[0];
      if (wf != 0)
      {
        if (domain->dimension >= 2)
        {
          s[1] = basis_function(r[1], ntype[in][1]);
          wf *= s[1];
        }
        else
          s[1] = 1;
        if (domain->dimension == 3 && wf != 0)
        {
          s[2] = basis_function(r[2], ntype[in][2]);
          wf *= s[2];
        }
        else
          s[2] = 1;
      }

      if (wf != 0)
      {
        if (solid.mat->rigid)
          solid.grid->rigid[in] = true;

        sd[0] = derivative_basis_function(r[0], ntype[in][0], inv_cellsize);
        if (domain->dimension >= 2)
          sd[1] = derivative_basis_function(r[1], ntype[in][1], inv_cellsize);
        if (domain->dimension == 3)
          sd[2] = derivative_basis_function(r[2], ntype[in][2], inv_cellsize);

        solid.neigh_n(ip, i) = in;
        solid.wf     (ip, i) = wf;

        if (domain->dimension == 3)
        {
          wfd[0] = sd[0] * s[1] * s[2];
          wfd[1] = s[0] * sd[1] * s[2];
          wfd[2] = s[0] * s[1] * sd[2];
        }
        else if (domain->dimension == 2)
        {
          wfd[0] = sd[0] * s[1];
          wfd[1] = s[0] * sd[1];
          wfd[2] = 0;
        }
        else
        {
          wfd[0] = sd[0];
          wfd[1] = 0;
          wfd[2] = 0;
        }

        solid.wfd(ip, i) = wfd;
      }
    }
  }

  if (update_Di && apic())
    solid.compute_inertia_tensor();

  //if (update->ntimestep == 0)
  //{
  //  // Reduce rigid_solids
  //  int rigid_solids_reduced = 0;

  //  MPI_Allreduce(&rigid_solids, &rigid_solids_reduced, 1, MPI_INT, MPI_LOR,
  //                universe->uworld);

  //  rigid_solids = rigid_solids_reduced;
  //}

  //if (rigid_solids)
  //  domain->grid->reduce_rigid_ghost_nodes();

  //update_Di = 0;
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
  grid.mass[in] = 0;
}

void Method::compute_mass_nodes(Solid &solid, int ip)
{
  for (int i = 0; i < solid.neigh_n.extent(1); i++)
  {
    int in = solid.neigh_n(ip, i);
    double wf = solid.wf(ip, i);

    if (!solid.grid->rigid[in] || solid.mat->rigid)
      solid.grid->mass[in] += wf*solid.mass[ip];
  }
}

void Method::reset_velocity_nodes(Grid &grid, int in)
{
  grid.v[in] = Vector3d();
  grid.mb[in] = Vector3d();
  if (temp)
    grid.T[in] = 0;
}

void Method::compute_velocity_nodes(Solid &solid, int ip)
{
  for (int i = 0; i < solid.neigh_n.extent(1); i++)
  {
    int in = solid.neigh_n(ip, i);
    double wf = solid.wf(ip, i);

    if (solid.grid->rigid[in] && !solid.mat->rigid)
      return;

    if (double node_mass = solid.grid->mass[in])
    {
      double normalized_wf = wf*solid.mass[ip]/node_mass;
      const Vector3d &dx = solid.grid->x0[in] - solid.x[ip];
    
      Vector3d vtemp = solid.v[ip];

      if (apic() || update->method->ge)
      {
        if (is_TL)
          vtemp += solid.Fdot[ip]*(solid.grid->x0[in] - solid.v[ip]);
        else
          vtemp += solid.L[ip]*(solid.grid->x0[in] - solid.x[ip]);
      }

      solid.grid->v[in] += normalized_wf*vtemp;

      if (solid.grid->rigid[in])
        solid.grid->mb[in] += normalized_wf*solid.v_update[ip];

      if (temp)
        solid.grid->T[in] += normalized_wf*solid.T[ip];
    }
  }
}

void Method::reset_force_nodes(Grid &grid, int in)
{
  grid.f[in] = Vector3d();
  grid.mb[in] = Vector3d();
  if (temp)
  {
    grid.Qint[in] = 0;
    grid.Qext[in] = 0;
  }
}

void Method::compute_force_nodes(Solid &solid, int ip)
{
  compute_internal_force_nodes(solid, ip);
  
  for (int i = 0; i < solid.neigh_n.extent(1); i++)
  {
    int in = solid.neigh_n(ip, i);
    double wf = solid.wf(ip, i);
    const Vector3d &wfd = solid.wfd(ip, i);

    if (!solid.grid->rigid[in])
      solid.grid->mb[in] += wf*solid.mbp[ip];

    if (temp)
    {
      if (solid.domain->axisymmetric)
      {
        error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
      }

      if (solid.grid->mass[in])
        solid.grid->Qext[in] += wf*solid.gamma[ip];

      solid.grid->Qint[in] += wfd.dot(solid.q[ip]);
    }
  }
}

void Method::update_grid_velocities(Grid &grid, int in)
{
  double T_update;

  Vector3d &v_update = grid.v_update[in] = grid.v[in];
  if (temp)
    T_update = grid.T_update[in] = grid.T[in];

  if (double mass = grid.mass[in])
  {
    if (!grid.rigid[in])
      v_update += update->dt*(grid.f[in] + grid.mb[in])/mass;

    if (temp)
      T_update += update->dt*(grid.Qint[in] + grid.Qext[in])/mass;
  }
}

void Method::compute_velocity_acceleration(Solid &solid, int ip)
{
  solid.v_update[ip] = Vector3d();
  solid.a[ip] = Vector3d();
  solid.f[ip] = Vector3d();

  for (int i = 0; i < solid.neigh_n.extent(1); i++)
  {
    int in = solid.neigh_n(ip, i);
    double wf = solid.wf(ip, i);

    solid.v_update[ip] += wf*solid.grid->v_update[in];

    if (solid.mat->rigid)
      continue;

    const Vector3d &delta_a = wf*(solid.grid->v_update[in] - solid.grid->v[in])/update->dt;
    solid.a[ip] += delta_a;
    solid.f[ip] += delta_a*solid.mass[ip];

    if (temp)
      solid.T[ip] += wf*solid.grid->T_update[in];
  }
}

void Method::update_position(Solid &solid, int ip)
{
  check_particle_in_domain(solid.x[ip] += update->dt*solid.v_update[ip], ip);
}

void Method::advance_particles(Solid &solid, int ip)
{
  Vector3d &v = solid.v[ip];

  if (update->sub_method_type == Update::SubMethodType::ASFLIP)
    solid.x[ip] += update->dt*v;

  v = (1 - update->PIC_FLIP)*solid.v_update[ip] + update->PIC_FLIP*(v + update->dt*solid.a[ip]);
}

void Method::compute_rate_deformation_gradient(bool doublemapping, Solid &solid, int ip)
{
  if (solid.mat->rigid)
    return;

  get_gradients(solid)[ip] = Matrix3d();
  if (temp)
    solid.q[ip] = Vector3d();
        
  Kokkos::View<Matrix3d*> &gradients = get_gradients(solid);
  const Kokkos::View<Vector3d*> &vn = doublemapping? solid.grid->v: solid.grid->v_update;
  
  for (int i = 0; i < solid.neigh_n.extent(1); i++)
  {
    int in = solid.neigh_n(ip, i);
    double wf = solid.wf(ip, i);
    const Vector3d &wfd = solid.wfd(ip, i);

    if (update->sub_method_type == Update::SubMethodType::APIC)
    {
      const Vector3d &dx = is_TL? solid.grid->x0[in] - solid.v[ip]:
                                  solid.grid->x[ip] - solid.grid->x0[in];

      Matrix3d gradient;
      for (int j = 0; j < domain->dimension; j++)
        for (int k = 0; k < domain->dimension; k++)
          gradient(j, k) += vn[in][j]*dx[k]*wf;

      gradients[ip] += gradient*solid.Di;
    }
    else
      for (int j = 0; j < domain->dimension; j++)
        for (int k = 0; k < domain->dimension; k++)
          gradients[ip](j, k) += vn[in][j]*wfd[k];

    if (domain->dimension == 2 && domain->axisymmetric)
      gradients[ip](2, 2) += vn[in][0]*wf/solid.v[ip][0];

    if (temp)
      solid.q[ip] -= wfd*(doublemapping? solid.grid->T: solid.grid->T_update)[in]
                        *(is_TL? solid.vol0: solid.vol)[ip]*solid.mat->invcp*solid.mat->kappa;
  }
}

void Method::update_deformation_gradient_determinant(Solid &solid, int ip)
{
  // FOR CPDI:
  //if (update->method->style == 1)
  //{
  //  solid.vol[ip] = 0.5*(solid.xpc[solid.nc*ip + 0][0]*solid.xpc[solid.nc*ip + 1][1] -
  //                          solid.xpc[solid.nc*ip + 1][0]*solid.xpc[solid.nc*ip + 0][1] +
  //                          solid.xpc[solid.nc*ip + 1][0]*solid.xpc[solid.nc*ip + 2][1] -
  //                          solid.xpc[solid.nc*ip + 2][0]*solid.xpc[solid.nc*ip + 1][1] +
  //                          solid.xpc[solid.nc*ip + 2][0]*solid.xpc[solid.nc*ip + 3][1] -
  //                          solid.xpc[solid.nc*ip + 3][0]*solid.xpc[solid.nc*ip + 2][1] +
  //                          solid.xpc[solid.nc*ip + 3][0]*solid.xpc[solid.nc*ip + 0][1] -
  //                          solid.xpc[solid.nc*ip + 0][0]*solid.xpc[solid.nc*ip + 3][1]);
  //  solid.J[ip] = solid.vol[ip]/solid.vol0[ip];
  //}
  //else
  //{

  solid.J[ip] = determinant(solid.F[ip]);
  solid.vol[ip] = solid.J[ip]*solid.vol0[ip];
}

// EB: remove this
#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333

void Method::update_deformation_gradient_stress(bool doublemapping, Solid &solid, int ip)
{
  if (solid.mat->rigid)
    return;

  update_deformation_gradient_matrix(solid, ip);

  solid.Finv[ip] = inverse(solid.F[ip]);

  update_deformation_gradient_determinant(solid, ip);

  if (solid.J[ip] <= 0.0 && solid.damage[ip] < 1.0)
  {
    cout << "Error: J[" << solid.ptag[ip] << "]<=0.0 == " << solid.J[ip] << endl;
    cout << "F[" << solid.ptag[ip] << "]:" << endl << solid.F[ip] << endl;
    cout << "Fdot[" << solid.ptag[ip] << "]:" << endl << solid.Fdot[ip] << endl;
    cout << "damage[" << solid.ptag[ip] << "]:" << endl << solid.damage[ip] << endl;
    error->one(FLERR, "");
  }

  solid.rho[ip] = solid.rho0[ip]/solid.J[ip];

  if (solid.mat->type == material->constitutive_model::LINEAR)
  {
    update_velocity_gradient_matrix(solid, ip);

    const Matrix3d &strain_increment = update->dt*solid.D[ip];
    solid.strain_el[ip] += strain_increment;
    solid.sigma[ip] += 2*solid.mat->G*strain_increment +
      solid.mat->lambda*strain_increment.trace()*Matrix3d::identity();

    if (is_TL)
      solid.vol0PK1[ip] = solid.vol0[ip]*solid.J[ip]*solid.R[ip]*solid.sigma[ip]*
                            solid.R[ip].transpose()*solid.Finv[ip].transpose();
    solid.gamma[ip] = 0;
  }
  else if (solid.mat->type == material->constitutive_model::NEO_HOOKEAN)
  {
    const Matrix3d &FinvT = solid.Finv[ip].transpose();
    const Matrix3d &PK1 = solid.mat->G*(solid.F[ip] - FinvT) + solid.mat->lambda*log(solid.J[ip])*FinvT;
    solid.vol0PK1[ip] = solid.vol0[ip]*PK1;
    solid.sigma[ip] = 1/solid.J[ip]*solid.F[ip]*PK1.transpose();

    solid.strain_el[ip] = 0.5*(solid.F[ip].transpose()*solid.F[ip] - Matrix3d::identity());
    solid.gamma[ip] = 0;
  }
  else
  {
    update_velocity_gradient_matrix(solid, ip);

    double pH = 0;
    double plastic_strain_increment = 0;
    Matrix3d sigma_dev;

    double T = solid.mat->cp? solid.T[ip]: 0;

    solid.mat->eos->compute_pressure(pH, solid.ienergy[ip], solid.J[ip], solid.rho[ip],
                                      solid.damage[ip], solid.D[ip], solid.grid->cellsize, T);

    if (solid.mat->cp)
      pH += solid.mat->alpha*(solid.T[ip] - solid.T0);

    sigma_dev = solid.mat->strength->update_deviatoric_stress(
      solid.sigma[ip], solid.D[ip], plastic_strain_increment,
      solid.eff_plastic_strain[ip], solid.eff_plastic_strain_rate[ip], solid.damage[ip],
      T);

    solid.eff_plastic_strain[ip] += plastic_strain_increment;

    // compute a characteristic time over which to average the plastic strain
    solid.eff_plastic_strain_rate[ip] += (plastic_strain_increment - solid.eff_plastic_strain_rate[ip]*update->dt)/
                                            1000/solid.grid->cellsize*solid.mat->signal_velocity;
    solid.eff_plastic_strain_rate[ip] = MAX(0.0, solid.eff_plastic_strain_rate[ip]);

    if (solid.mat->damage)
        solid.mat->damage->compute_damage(solid.damage_init[ip], solid.damage[ip], pH,
                                          sigma_dev, solid.eff_plastic_strain_rate[ip],
                                          plastic_strain_increment, T);

    if (solid.mat->temp)
    {
      solid.mat->temp->compute_heat_source(solid.T[ip], solid.gamma[ip], SQRT_3_OVER_2*sigma_dev.norm(),
                                           solid.eff_plastic_strain_rate[ip]);
      if (is_TL)
        solid.gamma[ip] *= solid.vol0[ip]*solid.mat->invcp;
      else
        solid.gamma[ip] *= solid.vol[ip]*solid.mat->invcp;
    }
    else
	  solid.gamma[ip] = 0;

    solid.sigma[ip] = -pH*(1 - (pH < 0? solid.damage[ip]: 0))*Matrix3d::identity() + sigma_dev;

      solid.strain_el[ip] =
        (update->dt*solid.D[ip].trace() + solid.strain_el[ip].trace())/3*Matrix3d::identity() +
        sigma_dev/solid.mat->G/(1 - (solid.damage[ip] > 1e-10? solid.damage[ip]: 0));

    if (is_TL)
      solid.vol0PK1[ip] = solid.vol0[ip]*solid.J[ip]*
        solid.R[ip]*solid.sigma[ip]*solid.R[ip].transpose()*
        solid.Finv[ip].transpose();
  }

  if (solid.damage[ip] >= 1.0)
    return;

  double p_wave_speed = sqrt((solid.mat->K + FOUR_THIRD*solid.mat->G)/solid.rho[ip]) +
                        MAX(MAX(fabs(solid.v[ip](0)), fabs(solid.v[ip](1))), fabs(solid.v[ip](2)));

  if (std::isnan(p_wave_speed))
  {
    cout << "Error: max_p_wave_speed is nan with ip=" << ip
      << ", ptag[ip]=" << solid.ptag[ip] << ", rho0[ip]=" << solid.rho0[ip]<< ", rho[ip]=" << solid.rho[ip]
      << ", K=" << solid.mat->K << ", G=" << solid.mat->G << ", J[ip]=" << solid.J[ip]
      << endl;
    error->one(FLERR, "");
  }
  else if (p_wave_speed < 0.0)
  {
    cout << "Error: p_wave_speed= " << p_wave_speed
      << " with ip=" << ip << ", rho[ip]=" << solid.rho[ip] << ", K=" << solid.mat->K
      << ", G=" << solid.mat->G << endl;
    error->one(FLERR, "");
  }
  
  double h_ratio = 1;

  if (is_TL)
  {
    Matrix3d eigenvalues = solid.F[ip];
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
        << "F=\n" <<  solid.F[ip] << endl
        << "eigenvalues of F:" << eigenvalues(0, 0) << "\t" << eigenvalues(1, 1) << "\t" << eigenvalues(2, 2) << endl;
      //cout << "esF.info()=" << esF.info() << endl;
      error->one(FLERR, "");
    }
  }
  else
  {
    h_ratio = 1;
  }

  solid.dtCFL[ip] = solid.grid->cellsize*h_ratio/p_wave_speed;

  if (std::isnan(solid.dtCFL[ip]))
  {
    cout << "Error: dtCFL = " << solid.dtCFL[ip] << "\n";
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
    for (int ip = 0; ip < domain->solids[isolid]->np_local; ip++)
    {
      dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL[ip]);
      if (!dtCFL)
      {
        cout << "Error: dtCFL == 0\n";
        cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
        error->one(FLERR, "");
      }
      else if (std::isnan(dtCFL))
      {
        cout << "Error: dtCFL = " << dtCFL << "\n";
        cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL[ip] << "\n";
        error->one(FLERR, "");
      }
    }

  MPI_Allreduce(&dtCFL, &dtCFL_reduced, 1, MPI_DOUBLE, MPI_MIN, universe->uworld);

  update->dt = dtCFL_reduced * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void Method::reset()
{
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
  {
    for (int ip = 0; ip < domain->solids[isolid]->np_local; ip++)
      domain->solids[isolid]->mbp[ip] = Vector3d();
  }
}
