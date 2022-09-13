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
#include <mpm_math.h>

using namespace std;

Method::Method(MPM *mpm):
  Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
  temp = false;
}

void Method::compute_grid_weight_functions_and_gradients(Solid &solid)
{
  bigint nnodes = solid.grid->nnodes;

  if (is_TL && update->atimestep || !nnodes)
    return;

  bool rigid = solid.mat->rigid;

  if (update->ntimestep == 0 && rigid)
    rigid_solids = 1;

  Kokkos::View<Vector3d*> &x0 = solid.grid->x0;
  Kokkos::View<Vector3i*> &ntype = solid.grid->ntype;
  Kokkos::View<tagint*> &map_ntag = solid.grid->map_ntag;
  Kokkos::View<bool*> grid_rigid = solid.grid->rigid;

  Kokkos::View<Vector3d*> x = solid.x;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  Kokkos::View<float**> wf = solid.wf;
  Kokkos::View<Vector3d**> wfd = solid.wfd;

  const float &inv_cellsize = 1/solid.grid->cellsize;
    
  int ny = solid.grid->ny_global;
  int nz = solid.grid->nz_global;

  Update::ShapeFunctions shape_function = update->shape_function;
  int half_support = shape_function == Update::ShapeFunctions::LINEAR? 1: 2;

  int dimension = domain->dimension;

  Vector3d boxlo(domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]);

  Kokkos::parallel_for("compute_grid_weight_functions_and_gradients", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    // Calculate what nodes particle ip will interact with:
    
    int count = 0;

    const Vector3d &xp = x[ip];
    const int &i0 = (xp[0] - boxlo[0])*inv_cellsize + 1 - half_support;
    const int &j0 = (xp[1] - boxlo[1])*inv_cellsize + 1 - half_support;
    const int &k0 = (xp[2] - boxlo[2])*inv_cellsize + 1 - half_support;

    for (int i = i0; i < i0 + 2*half_support; i++)
    {
      if (ny > 1)
      {
        for (int j = j0; j < j0 + 2*half_support; j++)
        {
          if (nz > 1)
          {
            for (int k = k0; k < k0 + 2*half_support; k++)
            {
              const int &tag = nz*ny*i + nz*j + k;

              if (tag < nnodes)
              {
                const tagint &in = map_ntag[tag];

                if (in != -1)
                  neigh_n(ip, count++) = in;
              }
            }
          }
          else
          {
            const int &tag = ny*i + j;

            if (tag < nnodes)
            {
              const tagint &in = map_ntag[tag];

              if (in != -1)
                neigh_n(ip, count++) = in;
            }
          }
        }
      }
      else if (i < nnodes)
        neigh_n(ip, count++) = i;
    }

    // not-flat-list grid version
    /*vector<vector<double>> s(dimension, vector<double>(2*half_support)), sd = s;

    for (int i = 0, power = 1; i < dimension; i++)
    {
      for (int j = 0; j < 2*half_support; j++)
      {
        const int &in = neigh_n(ip, j*power);
        
        const double &r = (xp - x0[in])[i]*inv_cellsize;
        
        if (shape_function == Update::ShapeFunctions::LINEAR)
        {
          s .at(i).at(j) = BasisFunction::           linear(r);
          sd.at(i).at(j) = BasisFunction::derivative_linear(r, inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::CUBIC_SPLINE)
        {
          s .at(i).at(j) = BasisFunction::           cubic_spline(r, ntype[in][j]);
          sd.at(i).at(j) = BasisFunction::derivative_cubic_spline(r, ntype[in][j], inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::QUADRATIC_SPLINE)
        {
          s .at(i).at(j) = BasisFunction::           quadratic_spline(r, ntype[in][j]);
          sd.at(i).at(j) = BasisFunction::derivative_quadratic_spline(r, ntype[in][j], inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::BERNSTEIN)
        {
          s .at(i).at(j) = BasisFunction::           bernstein_quadratic(r, ntype[in][j]);
          sd.at(i).at(j) = BasisFunction::derivative_bernstein_quadratic(r, ntype[in][j], inv_cellsize);
        }
      }

      power *= 2*half_support;
    }*/

    for (int i = 0; i < count; i++)
    {
      const int &in = neigh_n(ip, i);

      // Calculate the distance between each pair of particle/node:
      const Vector3d &r = (xp - x0[in])*inv_cellsize;

      Vector3d s, sd;

      for (int j = 0; j < 3; j++)
      {
        if (dimension <= j)
        {
          s [j] = 1;
          sd[j] = 0;
        }
        else if (shape_function == Update::ShapeFunctions::LINEAR)
        {
          s [j] = BasisFunction::           linear(r[j]);
          sd[j] = BasisFunction::derivative_linear(r[j], inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::CUBIC_SPLINE)
        {
          s [j] = BasisFunction::           cubic_spline(r[j], ntype[in][j]);
          sd[j] = BasisFunction::derivative_cubic_spline(r[j], ntype[in][j], inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::QUADRATIC_SPLINE)
        {
          s [j] = BasisFunction::           quadratic_spline(r[j], ntype[in][j]);
          sd[j] = BasisFunction::derivative_quadratic_spline(r[j], ntype[in][j], inv_cellsize);
        }
        else if (shape_function == Update::ShapeFunctions::BERNSTEIN)
        {
          s [j] = BasisFunction::           bernstein_quadratic(r[j], ntype[in][j]);
          sd[j] = BasisFunction::derivative_bernstein_quadratic(r[j], ntype[in][j], inv_cellsize);
        }
      }
      
      wf (ip, i)    = s [0]*s [1]*s [2];

      wfd(ip, i)[0] = sd[0]*s [1]*s [2];
      wfd(ip, i)[1] = s [0]*sd[1]*s [2];
      wfd(ip, i)[2] = s [0]*s [1]*sd[2];

      grid_rigid[in] = rigid;
    }
  });

  if (apic())
    solid.compute_inertia_tensor();

  if (update->ntimestep == 0)
  {
    // Reduce rigid_solids
    int rigid_solids_reduced = 0;

    MPI_Allreduce(&rigid_solids, &rigid_solids_reduced, 1, MPI_INT, MPI_LOR,
                  universe->uworld);

    rigid_solids = rigid_solids_reduced;
  }

  if (rigid_solids)
    domain->grid->reduce_rigid_ghost_nodes();

}

void Method::setup(vector<string> args)
{
  if (!args.empty())
    error->all(FLERR, "Illegal modify_method command: too many arguments.\n");

  if (update->sub_method_type == Update::SubMethodType::APIC || update->sub_method_type == Update::SubMethodType::MLS)
    update->PIC_FLIP = 0;
}

bool Method::apic()
{
  switch (update->sub_method_type)
  {
  case Update::SubMethodType::APIC:
  case Update::SubMethodType::MLS:
    return true;
  }

  return false;
}

vector<Grid *> Method::grids()
{
  vector<Grid *> grids;

  if (is_TL)
  {
    grids.reserve(domain->solids.size());

    for (Solid *solid: domain->solids)
      grids.push_back(solid->grid);
  }
  else
    grids.push_back(domain->grid);

  return grids;
}

void Method::reset_mass_nodes(Grid &grid)
{
  if (!is_TL || !update->atimestep)
    Kokkos::parallel_for("reset_mass_nodes", grid.nnodes_local + grid.nnodes_ghost,
    KOKKOS_LAMBDA (const int &in)
    {
      grid.mass[in] = 0;
    });
}

void Method::compute_mass_nodes(Solid &solid)
{
  Kokkos::View<float*> smass = solid.mass;
  Kokkos::View<float**> swf = solid.wf;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  bool srigid = solid.mat->rigid;

  Kokkos::View<float*> gmass = solid.grid->mass;
  Kokkos::View<bool*> grigid =  solid.grid->rigid;

  if (!is_TL || !update->atimestep)
    Kokkos::parallel_for("compute_mass_nodes", solid.neigh_policy,
    KOKKOS_LAMBDA (int ip, int i)
    {
      if (float wf = swf(ip, i))
      {
        int in = neigh_n(ip, i);

        if (!grigid[in] || srigid)
          Kokkos::atomic_add(&gmass[in], wf*smass[ip]);
      }
    });
}

void Method::reset_velocity_nodes(Grid &grid)
{
  bool temp = this->temp;
  Kokkos::View<Vector3d*> gv = grid.v;
  Kokkos::View<Vector3d*> gmb = grid.mb;
  Kokkos::View<float*> gT = grid.T;

  Kokkos::parallel_for("reset_velocity_nodes", grid.nnodes_local + grid.nnodes_ghost,
  KOKKOS_LAMBDA (const int &in)
  {
    gv[in] = Vector3d();
    gmb[in] = Vector3d();
    if (temp)
      gT[in] = 0;
  });
}

void Method::compute_velocity_nodes(Solid &solid)
{
  bool apic = this->apic();
  bool temp = this->temp;
  bool is_TL = this->is_TL;
  
  Kokkos::View<float*> smass = solid.mass;
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> sx0 = solid.x0;
  Kokkos::View<Vector3d*> sv = solid.v;
  Kokkos::View<Vector3d*> sv_update = solid.v_update;
  Kokkos::View<Matrix3d*> L = solid.L;
  Kokkos::View<Matrix3d*> Fdot = solid.Fdot;
  Kokkos::View<float*> sT = solid.T;
  Kokkos::View<float**> swf = solid.wf;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  bool srigid = solid.mat->rigid;
  
  Kokkos::View<float*> gmass = solid.grid->mass;
  Kokkos::View<bool*> grigid =  solid.grid->rigid;
  Kokkos::View<Vector3d*> gx0 = solid.grid->x0;
  Kokkos::View<Vector3d*> gv = solid.grid->v;
  Kokkos::View<Vector3d*> gmb = solid.grid->mb;
  Kokkos::View<float*> gT = solid.grid->T;

  Kokkos::parallel_for("compute_velocity_nodes", solid.neigh_policy,
  KOKKOS_LAMBDA (int ip, int i)
  {
    if (float wf = swf(ip, i))
    {
      int in = neigh_n(ip, i);

      if (grigid[in] && !srigid)
        return;

      if (float node_mass = gmass[in])
      {
        float normalized_wf = wf*smass[ip]/node_mass;
        const Vector3d &dx = gx0[in] - sx[ip];
    
        Vector3d vtemp = sv[ip];

        if (apic)
        {
          if (is_TL)
            vtemp += Fdot[ip]*(gx0[in] - sx0[ip]);
          else
            vtemp += L[ip]*(gx0[in] - sx[ip]);
        }

        gv[in].atomic_add(normalized_wf*vtemp);

        if (grigid[in])
          gmb[in].atomic_add(normalized_wf*sv_update[ip]);

        if (temp)
          Kokkos::atomic_add(&gT[in],  normalized_wf*sT[ip]);
      }
    }
  });
}

void Method::reset_force_nodes(Grid &grid)
{
  bool temp = this->temp;

  Kokkos::parallel_for("reset_force_nodes", grid.nnodes_local + grid.nnodes_ghost,
  KOKKOS_LAMBDA (const int &in)
  {
    grid.f[in] = Vector3d();
    grid.mb[in] = Vector3d();
    if (temp)
    {
      grid.Qint[in] = 0;
      grid.Qext[in] = 0;
    }
  });
}

void Method::compute_force_nodes(Solid &solid)
{
  bool temp = this->temp;
  bool is_TL = this->is_TL;
  Grid &grid = *solid.grid;
  Update::SubMethodType sub_method_type = update->sub_method_type;
  bool axisymmetric = domain->axisymmetric;

  Kokkos::View<float**> swf = solid.wf;
  Kokkos::View<Vector3d**> swfd = solid.wfd;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> sx0 = solid.x0;
  Kokkos::View<Matrix3d*> svol0PK1 = solid.vol0PK1;
  Kokkos::View<Vector3d*> smbp = solid.mbp;
  Kokkos::View<float*> sgamma = solid.gamma;
  Kokkos::View<Vector3d*> sq = solid.q;
  Kokkos::View<float*> svol = solid.vol;
  Kokkos::View<Matrix3d*> ssigma = solid.sigma;

  float &Di = solid.Di;

  Kokkos::View<float*> gmass = solid.grid->mass;
  Kokkos::View<Vector3d*> gx0 = solid.grid->x0;
  Kokkos::View<Vector3d*> gf = solid.grid->f;
  Kokkos::View<Vector3d*> gmb = solid.grid->mb;
  Kokkos::View<bool*> grigid =  solid.grid->rigid;
  Kokkos::View<float*> gQext = solid.grid->Qext;
  Kokkos::View<float*> gQint = solid.grid->Qint;

  Kokkos::parallel_for("compute_force_nodes", solid.neigh_policy,
  KOKKOS_LAMBDA (int ip, int i)
  {
    if (float wf = swf(ip, i))
    {
      int in = neigh_n(ip, i);
      if (gmass[in]) {
	const Vector3d &wfd = swfd(ip, i);
    
	Vector3d &f = gf[in];

	if (is_TL)
	  {
	    const Matrix3d &vol0PK1 = svol0PK1[ip];
	    const Vector3d &x0 = sx0[ip];

	    if (sub_method_type == Update::SubMethodType::MLS)
	      f.atomic_add(-vol0PK1*wf*Di*(grid.x0[in] - x0));
	    else
	      f.atomic_add(-vol0PK1*wfd);

	    if (axisymmetric)
	      Kokkos::atomic_add(&f[0], -vol0PK1(2, 2)*wf/x0[0]);
	  }
	else
	  {
	    const Matrix3d &vol_sigma = svol[ip]*ssigma[ip];
	    const Vector3d &x = sx[ip];

	    if (sub_method_type == Update::SubMethodType::MLS)
	      f.atomic_add(-vol_sigma*wf*Di*(gx0[in] - x));
	    else
	      f.atomic_add(-vol_sigma*wfd);

	    if (axisymmetric)
	      Kokkos::atomic_add(&f[0], -vol_sigma(2, 2)*wf/x[0]);
	  }

	if (!grigid[in])
	  gmb[in].atomic_add(wf*smbp[ip]);

	if (temp)
	  {
	    Kokkos::atomic_add(&gQint[in], wf*sgamma[ip] + wfd.dot(sq[ip]));
	  }
      }
    }
  });
}

void Method::update_grid_velocities(Grid &grid)
{
  float dt = update->dt;
  bool temp = this->temp;

  Kokkos::parallel_for("update_grid_velocities", grid.nnodes_local + grid.nnodes_ghost,
  KOKKOS_LAMBDA (const int &in)
  {
    float T_update;
  
    Vector3d &v_update = grid.v_update[in] = grid.v[in];
    if (temp)
      T_update = grid.T_update[in] = grid.T[in];

    if (float mass = grid.mass[in])
    {
      if (!grid.rigid[in])
        v_update += dt*(grid.f[in] + grid.mb[in])/mass;

      if (temp) {
        T_update += dt*(grid.Qint[in] + grid.Qext[in])/mass;
	grid.T_update[in] = T_update;
      }
    }
  });
}

void Method::compute_velocity_acceleration(Solid &solid)
{
  bool temp = this->temp;
  float dt = update->dt;
  Grid &grid = *solid.grid;
  bool rigid = solid.mat->rigid;

  Kokkos::View<float**> swf = solid.wf;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  Kokkos::View<Vector3d*> sv_update = solid.v_update;
  Kokkos::View<Vector3d*> sa = solid.a;
  Kokkos::View<Vector3d*> sf = solid.f;
  Kokkos::View<float*> sT = solid.T;
  Kokkos::View<float*> smass = solid.mass;
  
  Kokkos::parallel_for("compute_velocity_acceleration", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    Vector3d v_update;
    Vector3d a;
    Vector3d f;
    float T;
    if (temp)
      T = 0.9999*sT[ip];

    for (int i = 0; i < neigh_n.extent(1); i++)
    {
      if (float wf = swf(ip, i))
      {
        int in = neigh_n(ip, i);

        v_update += wf*grid.v_update[in];

        if (rigid)
          continue;

        const Vector3d &delta_a = wf*(grid.v_update[in] - grid.v[in])/dt;
        a += delta_a;
        f += delta_a*smass[ip];

        if (temp)
          T += wf*(grid.T_update[in] - grid.T[in]);
      }
    }

    sv_update[ip] = v_update;
    sa[ip] = a;
    sf[ip] = f;
    if (temp)
      sT[ip] = T;
  });
}

void Method::update_position(Solid &solid)
{
  float dt = update->dt;

  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> sv_update = solid.v_update;
  Kokkos::View<int*> serror_flag = solid.error_flag;

  Vector3d boxlo(domain->boxlo[0], domain->boxlo[1], domain->boxlo[2]);
  Vector3d boxhi(domain->boxhi[0], domain->boxhi[1], domain->boxhi[2]);

  int error_sum;
  if (!is_TL) {
    Kokkos::parallel_reduce("compute_velocity_acceleration0", solid.np_local,
			    KOKKOS_LAMBDA (const int &ip, int &error_)
    {
      sx[ip] += dt*sv_update[ip];

      if (sx[ip][0] < boxlo[0] ||
	  sx[ip][0] > boxhi[0] ||
	  sx[ip][1] < boxlo[1] ||
	  sx[ip][1] > boxhi[1] ||
	  sx[ip][2] < boxlo[2] ||
	  sx[ip][2] > boxhi[2])
	serror_flag[ip] = 1;
      error_ += serror_flag[ip];
    }, error_sum);

    if (error_sum) {
      Kokkos::View<Vector3d*>::HostMirror xp = create_mirror(solid.x);
      deep_copy(xp, solid.x);

      for (int ip = 0; ip < solid.np_local; ip++) {
        if (!domain->inside(xp[ip])) {
	  string msg = "Error: Particle " + to_string(ip) + " left the domain ("
	    + to_string(boxlo[0]) + "," + to_string(boxhi[0]) + ","
	    + to_string(boxlo[1]) + "," + to_string(boxhi[1]) + ","
	    + to_string(boxlo[2]) + "," + to_string(boxhi[2]) + "):\n"
	    + "x[" + to_string(ip) + "]=[" + to_string(xp[ip][0])
	    + ", " + to_string(xp[ip][1]) + ", "+ to_string(xp[ip][2]) + "]\n";

	  error->one(FLERR, msg);
	}
      }
    }
  } else {
    Kokkos::parallel_for("compute_velocity_acceleration0", solid.np_local,
			 KOKKOS_LAMBDA (const int &ip)
    {
      sx[ip] += dt*sv_update[ip];
    });
  }
}

void Method::advance_particles(Solid &solid)
{
  float PIC_FLIP = update->PIC_FLIP;
  float dt = update->dt;

  Kokkos::View<Vector3d*> sa = solid.a;
  Kokkos::View<Vector3d*> sv = solid.v;
  Kokkos::View<Vector3d*> sv_update = solid.v_update;

  Kokkos::parallel_for("advance_particles", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    Vector3d &v = sv[ip];

    v = (1 - PIC_FLIP)*sv_update[ip] + PIC_FLIP*(v + dt*sa[ip]);
  });
}

void Method::update_grid_positions(Grid &grid)
{
  if (is_TL)
  {
    float dt = update->dt;

    Kokkos::parallel_for("update_grid_positions", grid.nnodes_local + grid.nnodes_ghost,
    KOKKOS_LAMBDA (const int &in)
    {
      grid.x[in] += dt*grid.v[in];
    });
  }
}

void Method::compute_rate_deformation_gradient(bool doublemapping, Solid &solid)
{
  if (solid.mat->rigid)
    return;

  Kokkos::View<Matrix3d*> &gradients = is_TL? solid.Fdot: solid.L;
  const Kokkos::View<Vector3d*> &vn = doublemapping? solid.grid->v: solid.grid->v_update;

  bool temp = this->temp;
  Update::SubMethodType sub_method_type = update->sub_method_type;
  int dimension = domain->dimension;
  bool axisymmetric = domain->axisymmetric;
  bool is_TL = this->is_TL;
  Grid &grid = *solid.grid;
  float kappa_over_cp = solid.mat->kappa * solid.mat->invcp;

  Kokkos::View<float**> swf = solid.wf;
  Kokkos::View<Vector3d**> swfd = solid.wfd;
  Kokkos::View<int**> neigh_n = solid.neigh_n;
  Kokkos::View<Vector3d*> sx = is_TL? solid.x0: solid.x;
  Kokkos::View<Vector3d*> sv = solid.v;
  Kokkos::View<float*> svol = is_TL? solid.vol0: solid.vol;
  Kokkos::View<Vector3d*> sq = solid.q;
  Kokkos::View<float*> gT = doublemapping? grid.T: grid.T_update;
  float &Di = solid.Di;

  bool apic = sub_method_type == Update::SubMethodType::APIC? true: false;

  Kokkos::parallel_for("compute_rate_deformation_gradient0", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    Matrix3d gradient;
    Vector3d q;

    for (int i = 0; i < neigh_n.extent(1); i++)
    {
      if (float wf = swf(ip, i))
      {
        int in = neigh_n(ip, i);
        const Vector3d &wfd = apic? wf * (grid.x0[in] - sx [ip]) : swfd(ip, i);

	for (int j = 0; j < dimension; j++)
	  for (int k = 0; k < dimension; k++)
	    gradient(j, k) += vn[in][j]*wfd[k];

        if (axisymmetric)
          gradient(2, 2) += vn[in][0] * wf / sx[ip][0];

        if (temp)
          q -= swfd(ip, i)*gT[in]*svol[ip]*kappa_over_cp;
      }
    }

    gradients[ip] = apic? gradient*Di: gradient;

    if (temp)
      sq[ip] = q;
  });
}

// EB: remove this
#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333

void Method::update_deformation_gradient_stress(bool doublemapping, Solid &solid)
{
  if (solid.mat->rigid)
    return;

  bool is_TL = this->is_TL;
  float dt = update->dt;
  Grid &grid = *solid.grid;

  Kokkos::View<float*> svol = solid.vol;
  Kokkos::View<float*> svol0 = solid.vol0;
  Kokkos::View<float*> sJ = solid.J;
  Kokkos::View<Matrix3d*> sL = solid.L;
  Kokkos::View<Matrix3d*> sF = solid.F;
  Kokkos::View<Matrix3d*> sFdot = solid.Fdot;
  Kokkos::View<Matrix3d*> sFinv = solid.Finv;
  Kokkos::View<float*> srho = solid.rho;
  Kokkos::View<float*> srho0 = solid.rho0;

  Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    if (is_TL)
      sF[ip] += dt*sFdot[ip];
    else
      sF[ip] = (Matrix3d::identity() + dt*sL[ip])*sF[ip];

    sFinv[ip] = inverse(sF[ip]);
    
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

    sJ[ip] = determinant(sF[ip]);
    svol[ip] = sJ[ip]*svol0[ip];

    /*if (solid.J[ip] <= 0.0 && solid.damage[ip] < 1.0)
    {
      cout << "Error: J[" << solid.ptag[ip] << "]<=0.0 == " << solid.J[ip] << endl;
      cout << "F[" << solid.ptag[ip] << "]:" << endl << solid.F[ip] << endl;
      cout << "Fdot[" << solid.ptag[ip] << "]:" << endl << solid.Fdot[ip] << endl;
      cout << "damage[" << solid.ptag[ip] << "]:" << endl << solid.damage[ip] << endl;
      error->one(FLERR, "");
    }*/

    srho[ip] = srho0[ip]/sJ[ip];
  });

  if (solid.mat->type != Material::constitutive_model::NEO_HOOKEAN)
  {
    Kokkos::View<Matrix3d*> sD = solid.D;
    Kokkos::View<Matrix3d*> sR = solid.R;

    Kokkos::parallel_for("update_deformation_gradient_stress1", solid.np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      if (is_TL)
      {
        // polar decomposition of the deformation gradient, F = R*U
        if (!MPM_Math::PolDec(sF[ip], sR[ip]))
        {
          /*cout << "Polar decomposition of deformation gradient failed for particle " << ip << ".\n";
          cout << "F:" << endl << sF[ip] << endl;
          cout << "timestep:" << endl << update->ntimestep << endl;
          error->one(FLERR, "");*/
        }

        // In TLMPM. L is computed from Fdot:
        sL[ip] = sFdot[ip]*sFinv[ip];
        sD[ip] = 0.5*(sR[ip].transpose()*(sL[ip] + sL[ip].transpose())*sR[ip]);
      }
      else
        sD[ip] = 0.5*(sL[ip] + sL[ip].transpose());
    });
  }

  float G = solid.mat->G;
  float lambda = solid.mat->lambda;

  if (solid.mat->type == Material::constitutive_model::LINEAR)
  {
    Kokkos::View<Matrix3d*> svol0PK1 = solid.vol0PK1;
    Kokkos::View<Matrix3d*> sR = solid.R;
    Kokkos::View<Matrix3d*> sD = solid.D;
    Kokkos::View<Matrix3d*> sstrain_el = solid.strain_el;
    Kokkos::View<Matrix3d*> ssigma = solid.sigma;
    Kokkos::View<float*> sgamma = solid.gamma;

    Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      const Matrix3d &strain_increment = dt*sD[ip];
      sstrain_el[ip] += strain_increment;
      ssigma[ip] += 2*G*strain_increment +
        lambda*strain_increment.trace()*Matrix3d::identity();

      if (is_TL)
        svol0PK1[ip] = svol0[ip]*sJ[ip]*sR[ip]*ssigma[ip]*
	  sR[ip].transpose()*sFinv[ip].transpose();
      sgamma[ip] = 0;
    });
  }
  else if (solid.mat->type == Material::constitutive_model::NEO_HOOKEAN)
  {
    Kokkos::View<Matrix3d*> svol0PK1 = solid.vol0PK1;
    Kokkos::View<Matrix3d*> sstrain_el = solid.strain_el;
    Kokkos::View<Matrix3d*> ssigma = solid.sigma;
    Kokkos::View<float*> sgamma = solid.gamma;
    Kokkos::View<float*> sT = solid.T;
    float alpha = solid.mat->alpha;
    float K = solid.mat->K;
    float T0 = solid.T0;

    Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      const Matrix3d &FinvT = sFinv[ip].transpose();
      const Matrix3d &PK1 = G*(sF[ip] - FinvT) + lambda*Kokkos::Experimental::log(sJ[ip])*FinvT;
      svol0PK1[ip] = svol0[ip]*PK1;
      ssigma[ip] = 1/sJ[ip]*sF[ip]*PK1.transpose();

      sstrain_el[ip] = 0.5*(sF[ip].transpose()*sF[ip] - Matrix3d::identity());
      sgamma[ip] = 0;
      if (alpha != 0) {
	ssigma[ip] -= 3 * K * alpha * (sT[ip] - T0) * Matrix3d::identity();
	sstrain_el[ip] -= alpha * (sT[ip] - T0) * Matrix3d::identity();
      }
    });
  }
  else
  {
    Kokkos::View<float*> pH("pH", solid.np_local);
    Kokkos::View<float*> plastic_strain_increment("plastic_strain_increment", solid.np_local);
    Kokkos::View<Matrix3d*> sigma_dev("sigma_dev", solid.np_local);

    solid.mat->eos->compute_pressure(solid, pH);

    float alpha = solid.mat->alpha;
    Kokkos::View<float*> sT = solid.T;
    float &sT0 = solid.T0;

    if (solid.mat->cp)
      Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
      KOKKOS_LAMBDA (const int &ip)
      {
        pH[ip] += alpha*(sT[ip] - sT0);
      });
      
    solid.mat->strength->update_deviatoric_stress(solid, plastic_strain_increment, sigma_dev);
      
    float signal_velocity = solid.mat->signal_velocity;
    Kokkos::View<float*> seff_plastic_strain = solid.eff_plastic_strain;
    Kokkos::View<float*> seff_plastic_strain_rate = solid.eff_plastic_strain_rate;

    Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      seff_plastic_strain[ip] += plastic_strain_increment[ip];

      // compute a characteristic time over which to average the plastic strain
      seff_plastic_strain_rate[ip] += (plastic_strain_increment[ip] - seff_plastic_strain_rate[ip]*dt)/
                                              1000/grid.cellsize*signal_velocity;
      seff_plastic_strain_rate[ip] = MAX(0.0, seff_plastic_strain_rate[ip]);
    });

    if (solid.mat->damage)
      solid.mat->damage->compute_damage(solid, pH, sigma_dev, plastic_strain_increment);


    float invcp = solid.mat->invcp;
    Kokkos::View<float*> sgamma = solid.gamma;

    if (solid.mat->temp)
    {
      solid.mat->temp->compute_heat_source(solid, sigma_dev);
      
      Kokkos::View<float*> svol = solid.vol;
      Kokkos::View<float*> svol0 = solid.vol0;

      Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
      KOKKOS_LAMBDA (const int &ip)
      {
        if (is_TL)
          sgamma[ip] *= svol0[ip]*invcp;
        else
          sgamma[ip] *= svol[ip]*invcp;
      });
    }
    else
      Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
      KOKKOS_LAMBDA (const int &ip)
      {
	    sgamma[ip] = 0;
      });
    
    Kokkos::View<Matrix3d*> sstrain_el = solid.strain_el;
    Kokkos::View<Matrix3d*> ssigma = solid.sigma;
    Kokkos::View<float*> sdamage = solid.damage;
    Kokkos::View<Matrix3d*> sR = solid.R;
    Kokkos::View<Matrix3d*> sD = solid.D;
    Kokkos::View<Matrix3d*> svol0PK1 = solid.vol0PK1;

    Kokkos::parallel_for("update_deformation_gradient_stress0", solid.np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      ssigma[ip] = -pH[ip]*(1 - (pH[ip] < 0? sdamage[ip]: 0))*Matrix3d::identity() + sigma_dev[ip];

        sstrain_el[ip] =
          (dt*sD[ip].trace() + sstrain_el[ip].trace())/3*Matrix3d::identity() +
          sigma_dev[ip]/G/(1 - (sdamage[ip] > 1e-10? sdamage[ip]: 0));

      if (is_TL)
        svol0PK1[ip] = svol0[ip]*sJ[ip]*
          sR[ip]*ssigma[ip]*sR[ip].transpose()*
          sFinv[ip].transpose();
    });
  }

  float K = solid.mat->K;

  Kokkos::View<Vector3d*> sv = solid.v;
  Kokkos::View<float*> sdtCFL = solid.dtCFL;
  Kokkos::View<float*> sdamage = solid.damage;

  Kokkos::parallel_for("update_deformation_gradient_stress2", solid.np_local,
  KOKKOS_LAMBDA (const int &ip)
  {
    if (sdamage[ip] >= 1)
      return;
      
    float p_wave_speed = Kokkos::Experimental::sqrt((K + G*4/3)/srho[ip]) +
                          MAX(MAX(Kokkos::Experimental::abs(sv[ip](0)),
                                  Kokkos::Experimental::abs(sv[ip](1))),
                                  Kokkos::Experimental::abs(sv[ip](2)));

    /*if (std::isnan(p_wave_speed))
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
    }*/
  
    float h_ratio = 1;

    if (is_TL)
    {
      Matrix3d eigenvalues = sF[ip];
      eigendecompose(eigenvalues);
      // EB: revisit
      /*if (esF.info()== Successfalse)
      {
        h_ratio = MIN(h_ratio, fabs(eigenvalues(0, 0)));
        h_ratio = MIN(h_ratio, fabs(eigenvalues(1, 1)));
        h_ratio = MIN(h_ratio, fabs(eigenvalues(2, 2)));
      }*/

      /*if (h_ratio == 0)
      {
        cout << "min_h_ratio == 0 with ip=" << ip
          << "F=\n" <<  solid.F[ip] << endl
          << "eigenvalues of F:" << eigenvalues(0, 0) << "\t" << eigenvalues(1, 1) << "\t" << eigenvalues(2, 2) << endl;
        //cout << "esF.info()=" << esF.info() << endl;
        error->one(FLERR, "");
      }*/
    }
    else
    {
      h_ratio = 1;
    }

    sdtCFL[ip] = grid.cellsize*h_ratio/p_wave_speed;

    /*if (std::isnan(solid.dtCFL[ip]))
    {
      cout << "Error: dtCFL = " << solid.dtCFL[ip] << "\n";
      cout << "p_wave_speed = " << p_wave_speed
        << ", grid->cellsize=" << solid.grid->cellsize << endl;
      error->one(FLERR, "");
    }*/
  });
}

void Method::adjust_dt()
{
  if (update->dt_constant)
    return; // dt is set as a constant, do not update

  float dtCFL = 1.0e22;
  float dtCFL_reduced = 1.0e22;

  for (Solid *solid: domain->solids)
  {
    Kokkos::View<float*> solid_dtCFL = solid->dtCFL;

    Kokkos::parallel_reduce("update_deformation_gradient_stress2", solid->np_local,
    KOKKOS_LAMBDA (const int &ip, float &dtCFL1)
    {
      dtCFL1 = MIN(dtCFL1, solid_dtCFL[ip]);
      /*if (!dtCFL)
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
      }*/
    }, Kokkos::Min<float>(dtCFL));
  }

  MPI_Allreduce(&dtCFL, &dtCFL_reduced, 1, MPI_FLOAT, MPI_MIN, universe->uworld);

  update->dt = dtCFL_reduced*update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void Method::reset()
{
  for (Solid *solid: domain->solids)
  {
    Kokkos::View<Vector3d*> mbp = solid->mbp;
    
    Kokkos::parallel_for("reset", solid->np_local,
    KOKKOS_LAMBDA (const int &ip)
    {
      mbp[ip] = Vector3d();
    });
  }
}

void Method::exchange_particles()
{
  if (is_TL || universe->nprocs < 2)
    return;
  
  /*
    HOW TO PARALLELIZE THIS:
    kokkos view not_in_domain_indexes

    parallel_scan i in all_particles:
      if (not_in_domain)
        cumulative_not_in_domain++;
      if (fimal)
        not_in_domain_indexes[cumulative_not_in_domain] = i

    parallel_for i in not_in_domain_indexes
      index0 = not_in_domain_indexes[i]
      index1 = length(all_particles) - (length(not_in_domain_indexes) - i - 1)

      if indexes not equal
        swap particles at index0 and index1

    pack and send last length(not_in_domain_indexes) particles, and shorten all_particles
  */

  int ip;
  // vector<int> np_send;
  vector<vector<float>> buf_send_vect(universe->nprocs);
  vector<vector<float>> buf_recv_vect(universe->nprocs);
  vector<int> unpack_list;
  int owner = 0;

  // Identify the particles that are not in the subdomain
  // and transfer their variables to the buffer:

  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    for (int iproc = 0; iproc < universe->nprocs; iproc++) {
      buf_send_vect[iproc].clear();
      buf_recv_vect[iproc].clear();
    }
    
    Kokkos::View<Vector3d*>::HostMirror xp = create_mirror(domain->solids[isolid]->x);
    deep_copy(xp, domain->solids[isolid]->x);

    // np_send.assign(universe->nprocs, 0);

    ip = 0;
    while (ip < domain->solids[isolid]->np_local) {
      owner = domain->which_CPU_owns_me(xp[ip][0], xp[ip][1], xp[ip][2]);
      if (owner != universe->me) {
        // The particle is not located in the subdomain anymore:
        // transfer it to the buffer
        domain->solids[isolid]->pack_particle(ip, buf_send_vect[owner]);
        domain->solids[isolid]->copy_particle(
            domain->solids[isolid]->np_local - 1, ip);
        domain->solids[isolid]->np_local--;
	// np_send[owner]++;
      } else {
        ip++;
      }
    }

    // if (np_local_old != domain->solids[isolid]->np_local) {
    //   domain->solids[isolid]->grow(domain->solids[isolid]->np_local);
    // }

    // New send and receive
    int size_r, size_s;
    int jproc;
    for (int i = 0; i < universe->sendnrecv.size(); i++) {
      if (universe->sendnrecv[i][0] == 0) {
	// Receive
	jproc = universe->sendnrecv[i][1];

        MPI_Recv(&size_r, 1, MPI_INT, jproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

	if (size_r > 0) {
	  buf_recv_vect[jproc].resize(size_r);
          MPI_Recv(&buf_recv_vect[jproc][0], size_r, MPI_FLOAT, jproc, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      } else {
	// Send
	jproc = universe->sendnrecv[i][1];
	size_s = buf_send_vect[jproc].size();
        MPI_Send(&size_s, 1, MPI_INT, jproc, 0, MPI_COMM_WORLD);

	if (size_s > 0) {
          MPI_Send(buf_send_vect[jproc].data(), size_s, MPI_FLOAT, jproc, 0,
                   MPI_COMM_WORLD);
        }
      }
    }


    // Check what particles are within the subdomain:
    for (int sproc = 0; sproc < universe->nprocs; sproc++) {
      if (sproc != universe->me) {
        unpack_list.clear();
        ip = 0;
        while (ip < buf_recv_vect[sproc].size()) {
          if (domain->inside_subdomain(buf_recv_vect[sproc][ip + 1],
                                       buf_recv_vect[sproc][ip + 2],
                                       buf_recv_vect[sproc][ip + 3])) {
            unpack_list.push_back(ip);
          } else {
            error->one(FLERR, "Particle received from CPU" + to_string(sproc) +
                                  " that did not belong in this CPU.\n ");
          }
          ip += domain->solids[isolid]->comm_n;
        }

        if (unpack_list.size() > 0) {
          domain->solids[isolid]->grow(domain->solids[isolid]->np_local +
                                       unpack_list.size());

          // Unpack buffer:
          domain->solids[isolid]->unpack_particle(
              domain->solids[isolid]->np_local, unpack_list,
              buf_recv_vect[sproc]);
        }
      }
    }
  }
}
