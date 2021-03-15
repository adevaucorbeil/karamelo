/* ----------------------------------------------------------------------
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

#include "solid.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "material.h"
#include "memory.h"
#include "method.h"
#include "mpm.h"
#include "mpm_math.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues> 
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <string>
#include <vector>

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

using namespace std;
using namespace Eigen;
using namespace MPM_Math;


#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333

Solid::Solid(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == NULL)
  {
    error->all(FLERR, "Error: a method should be defined before creating a solid!\n");
  }

  if (args.size() < 2)
  {
    string error_str = "Error: solid command not enough arguments. ";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (usage.find(args[1]) == usage.end())
  {
    string error_str = "Error, keyword \033[1;31m" + args[1] + "\033[0m unknown!\n";
    for (auto &x : usage)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (args.size() < Nargs.find(args[1])->second)
  {
    error->all(FLERR, "Error: not enough arguments.\n"
	       + usage.find(args[1])->second);
  }

  cout << "Creating new solid with ID: " << args[0] << endl;

  method_type = update->method_type;
  id          = args[0];

  np = 0;

  if (update->method->is_CPDI) {
    nc = pow(2, domain->dimension);
  }
  else
    nc = 0;

  mat = NULL;

  if (update->method->is_TL) {
    is_TL = true;
    grid = new Grid(mpm);
  }
  else {
    is_TL = false;
    grid = domain->grid;
  }

  if (update->method->method_type.compare("APIC") == 0) {
    apic = true;
  } else {
    apic = false;
  }

  dtCFL = 1.0e22;
  vtot  = 0;
  mtot = 0;

  if (args[1].compare("region") == 0)
  {
    // Set material, cellsize, and initial temperature:
    options(&args, args.begin() + 4);

    // Create particles:
    populate(args);
  }
  else if (args[1].compare("mesh") == 0)
  {
    // Set material and cellsize and initial temperature:
    options(&args, args.begin() + 3);

    read_mesh(args[2]);
  }
  else if (args[1].compare("file") == 0)
  {
    // Set material and cellsize and initial temperature:
    options(&args, args.begin() + 3);

    read_file(args[2]);
  }

  comm_n = 54; // Number of double to pack for particle exchange between CPUs.

}

Solid::~Solid()
{
  if (is_TL) delete grid;
}

void Solid::init()
{
  cout << "Bounds for " << id << ":\n";
  cout << "xlo xhi: " << solidlo[0] << " " << solidhi[0] << endl;
  cout << "ylo yhi: " << solidlo[1] << " " << solidhi[1] << endl;
  cout << "zlo zhi: " << solidlo[2] << " " << solidhi[2] << endl;

  // Calculate total volume:
  double vtot_local = 0;
  double mtot_local = 0;
  for (int ip=0; ip<np_local; ip++)
    {
      vtot_local += vol[ip];
      mtot_local += mass[ip];
  }

  MPI_Allreduce(&vtot_local, &vtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&mtot_local, &mtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

  cout << "Solid " << id << " total volume = " << vtot << endl;
  cout << "Solid " << id << " total mass = " << mtot << endl;

  if (grid->nnodes == 0) grid->init(solidlo, solidhi);

  if (np == 0) {
    error->one(FLERR,"Error: solid does not have any particles.\n");
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
  if (args->end() < it + 3)
  {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it)
  {
    int iMat = material->find_material(*it);

    if (iMat == -1)
    {
      cout << "Error: could not find material named " << *it << endl;
      error->all(FLERR,"\n");
    }

    mat = &material->materials[iMat]; // point mat to the right material

    it++;

    grid->setup(*it); // set the grid cellsize

    it++;
    T0 = input->parsev(*it); // set initial temperature

    it++;

    if (it != args->end())
    {
      string error_str = "Error: too many arguments\n";
      for (auto &x : usage)
        error_str += x.second;
    }
  }
}

void Solid::grow(int nparticles)
{
  ptag.resize(nparticles);
  x0.resize(nparticles);
  x.resize(nparticles);

  if (method_type.compare("tlcpdi") == 0
      || method_type.compare("ulcpdi") == 0)
    {

      if (update->method->style == 0)
	{ // CPDI-R4
	  rp0.resize(domain->dimension*nparticles);
	  rp.resize(domain->dimension*nparticles);
	}
      if (update->method->style == 1)
	{ // CPDI-Q4
	  xpc0.resize(nc * nparticles);
	  xpc.resize(nc * nparticles);
	}
    }

  if (method_type.compare("tlcpdi2") == 0
      || method_type.compare("ulcpdi2") == 0)
    {
      xpc0.resize(nparticles);
      xpc.resize(nparticles);
    }

  v.resize(nparticles);
  v_update.resize(nparticles);
  a.resize(nparticles);
  mbp.resize(nparticles);
  f.resize(nparticles);
  sigma.resize(nparticles);
  strain_el.resize(nparticles);
  vol0PK1.resize(nparticles);
  L.resize(nparticles);
  F.resize(nparticles);
  R.resize(nparticles);
  D.resize(nparticles);
  Finv.resize(nparticles);
  Fdot.resize(nparticles);
  if (apic)
    Di.resize(nparticles);

  vol0.resize(nparticles);
  vol.resize(nparticles);
  rho0.resize(nparticles);
  rho.resize(nparticles);
  mass.resize(nparticles);
  eff_plastic_strain.resize(nparticles);
  eff_plastic_strain_rate.resize(nparticles);
  damage.resize(nparticles);
  damage_init.resize(nparticles);
  ienergy.resize(nparticles);
  mask.resize(nparticles);
  J.resize(nparticles);

  numneigh_pn.resize(nparticles);
  neigh_pn.resize(nparticles);
  wf_pn.resize(nparticles);
  if (nc != 0)
    wf_pn_corners.resize(nc * nparticles);
  wfd_pn.resize(nparticles);

  bigint nnodes = grid->nnodes_local + grid->nnodes_ghost;

  numneigh_np.resize(nnodes);
  neigh_np.resize(nnodes);
  wf_np.resize(nnodes);
  wfd_np.resize(nnodes);

  if (update->method->temp) {
    T.resize(nparticles);
    gamma.resize(nparticles);
    q.resize(nparticles);
  }
}

void Solid::compute_mass_nodes(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
    {
      if (reset) grid->mass[in] = 0;

      if (grid->rigid[in] && !mat->rigid) continue;

      for (int j = 0; j < numneigh_np[in]; j++)
	{
	  ip = neigh_np[in][j];
	  grid->mass[in] += wf_np[in][j] * mass[ip];
	}
    }
  return;
}

void Solid::compute_velocity_nodes(bool reset)
{
  Eigen::Vector3d vtemp, vtemp_update;
  //double mass_rigid;
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
    {
      grid->v[in].setZero();
      //grid->v_update[in].setZero();
      if (grid->rigid[in]) {
	grid->mb[in].setZero();
      }
    }

    if (grid->rigid[in] && !mat->rigid) continue;

    if (grid->mass[in] > 0)
    {
      vtemp.setZero();
      if (grid->rigid[in])
	vtemp_update.setZero();

      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
	if (grid->rigid[in]) {
	  vtemp_update += (wf_np[in][j] * mass[ip]) * v_update[ip];
	}
	if (update->method->ge) {
	  vtemp += (wf_np[in][j] * mass[ip]) *
	    (v[ip] + L[ip] * (grid->x0[in] - x[ip]));
	} else {
	  vtemp += wf_np[in][j] * mass[ip] * v[ip];
	}
        // grid->v[in] += (wf_np[in][j] * mass[ip]) * v[ip]/ grid->mass[in];
      }
      vtemp /= grid->mass[in];
      grid->v[in] += vtemp;
      if (grid->rigid[in]) {
	vtemp_update /= grid->mass[in];
	grid->mb[in] += vtemp_update; // This should be grid->v_update[in], but we are using mb to make the reduction of ghost particles easy. It will be copied to grid->v_update[in] in Grid::update_grid_velocities()
      }
      // if (isnan(grid->v_update[in][0]))
      //   cout << "in=" << in << "\tvn=[" << grid->v[in][0] << ", " << grid->v[in][1]
      //        << ", " << grid->v[in][2] << "]\tvp=[" << v[ip][0] << ", " << v[ip][1]
      //        << ", " << v[ip][2] << "],\tvn_update=[" << grid->v_update[in][0]
      //        << ", " << grid->v_update[in][1] << ", " << grid->v_update[in][2] << "]\n";
    }
  }
}

void Solid::compute_velocity_nodes_APIC(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
      grid->v[in].setZero();

    if (grid->rigid[in] && !mat->rigid)
      continue;

    if (grid->mass[in] > 0)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        grid->v[in] += (wf_np[in][j] * mass[ip]) *
	  (v[ip] + Fdot[ip] * (grid->x0[in] - x0[ip]));
      }
      grid->v[in] /= grid->mass[in];
    }
  }
}

void Solid::compute_external_forces_nodes(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
      grid->mb[in].setZero();

    if (grid->rigid[in])
      continue;

    if (grid->mass[in] > 0)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        grid->mb[in] += wf_np[in][j] * mbp[ip];
      }
    }
  }
}

void Solid::compute_internal_forces_nodes_TL()
{
  Eigen::Vector3d ftemp;
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (grid->rigid[in])
    {
      grid->f[in].setZero();
      continue;
    }

    ftemp.setZero();
    for (int j = 0; j < numneigh_np[in]; j++)
    {
      ip = neigh_np[in][j];
      ftemp -= vol0PK1[ip] * wfd_np[in][j];

      if (domain->axisymmetric == true)
      {
        ftemp[0] -= vol0PK1[ip](2, 2) * wf_np[in][j] / x0[ip][0];
      }
    }

    grid->f[in] = ftemp;
  }
}

void Solid::compute_internal_forces_nodes_UL(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
    {
      if (reset) grid->f[in].setZero();
      for (int j = 0; j < numneigh_np[in]; j++)
	{
	  ip = neigh_np[in][j];
	  grid->f[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
	}

      if (domain->axisymmetric == true)
	{
	  for (int j = 0; j < numneigh_np[in]; j++)
	    {
	      ip = neigh_np[in][j];
	      grid->f[in][0] -= vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
	    }
	}
    }
}

void Solid::compute_particle_velocities_and_positions()
{

  vector<Eigen::Vector3d> vc_update;
  vc_update.resize(nc);

  int in;

  bool update_corners;

  if ((method_type.compare("tlcpdi") == 0 ||
       method_type.compare("ulcpdi") == 0) &&
      (update->method->style == 1))
  {
    update_corners = true;
  }
  else
    update_corners = false;

  for (int ip = 0; ip < np_local; ip++)
  {
    v_update[ip].setZero();
    if (update_corners)
      for (int i = 0; i < nc; i++)
        vc_update[i].setZero();

    for (int j = 0; j < numneigh_pn[ip]; j++)
      {
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * grid->v_update[in];
      x[ip] += update->dt * wf_pn[ip][j] * grid->v_update[in];
      // if (isnan(x[ip](0)))
      //   cout << "ip=" << ip << "\tx=[" << x[ip](0) << "," << x[ip](1) << ","
      //        << x[ip](2) << "]\tin=" << in << "\tvn_update=["
      //        << grid->v_update[in](0) << "," << grid->v_update[in](1) << ","
      //        << grid->v_update[in](2) << "]\twf_pn=" << wf_pn[ip][j] << "\n";

      if (update_corners)
      {
        for (int ic = 0; ic < nc; ic++)
        {
          vc_update[ic] += wf_pn_corners[nc * ip + ic][j] * grid->v_update[in];
        }
      }
    }

    if (!is_TL)
    {
      // Check if the particle is within the box's domain:
      if (domain->inside(x[ip]) == 0)
      {
        cout << "Error: Particle " << ip << " left the domain ("
             << domain->boxlo[0] << "," << domain->boxhi[0] << ","
             << domain->boxlo[1] << "," << domain->boxhi[1] << ","
             << domain->boxlo[2] << "," << domain->boxhi[2] << ",):\n"
             << x[ip] << endl;
        error->one(FLERR, "");
      }
    }

    if (update_corners)
    {
      for (int ic = 0; ic < nc; ic++)
      {
        xpc[nc * ip + ic] += update->dt * vc_update[ic];
      }
    }
  }
}

void Solid::compute_particle_acceleration()
{
  double inv_dt = 1.0/update->dt;

  int in;

  for (int ip = 0; ip < np_local; ip++){
    a[ip].setZero();
    if (mat->rigid)
      continue;
    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      a[ip] += wf_pn[ip][j] * (grid->v_update[in] - grid->v[in]);
      // if (ptag[ip] == 101) {
      // 	printf("dv[%d]=[%4.3e %4.3e %4.3e], ", grid->ntag[in], (grid->v_update[in][0] - grid->v[in][0])*inv_dt, (grid->v_update[in][1] - grid->v[in][1])*inv_dt, (grid->v_update[in][2] - grid->v[in][2])*inv_dt);
      // 	cout << "rigid=" << grid->rigid[in] << endl;
      // }
    }
    a[ip] *= inv_dt;
    f[ip] = a[ip] * mass[ip];
  }
}


void Solid::update_particle_velocities(double alpha)
{
  for (int ip = 0; ip < np_local; ip++)
    {
      v[ip] = (1 - alpha) * v_update[ip] + alpha * (v[ip] + update->dt * a[ip]);
      // if (ptag[ip] == 101) {
      // 	printf("v=[%4.3e %4.3e %4.3e]\tv_update=[%4.3e %4.3e %4.3e]\ta=[%4.3e %4.3e %4.3e]\n", v[ip][0], v[ip][1], v[ip][2], v_update[ip][0], v_update[ip][1], v_update[ip][2], a[ip][0], a[ip][1], a[ip][2]);
      // }
  }
}

void Solid::compute_rate_deformation_gradient_TL() {
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *vn = &grid->v;

  if (domain->dimension == 1) {
    for (int ip = 0; ip < np_local; ip++) {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
      }
    }
  } else if ((domain->dimension == 2) && (domain->axisymmetric == true)) {
    for (int ip = 0; ip < np_local; ip++) {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
        Fdot[ip](2, 2) += (*vn)[in][0] * wf_pn[ip][j] / x0[ip][0];
      }
    }
  } else if ((domain->dimension == 2) && (domain->axisymmetric == false)) {
    for (int ip = 0; ip < np_local; ip++) {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip = 0; ip < np_local; ip++) {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](0, 2) += (*vn)[in][0] * wfd_pn[ip][j][2];
        Fdot[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
        Fdot[ip](1, 2) += (*vn)[in][1] * wfd_pn[ip][j][2];
        Fdot[ip](2, 0) += (*vn)[in][2] * wfd_pn[ip][j][0];
        Fdot[ip](2, 1) += (*vn)[in][2] * wfd_pn[ip][j][1];
        Fdot[ip](2, 2) += (*vn)[in][2] * wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_MUSL()
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *vn = &grid->v;

  if (domain->dimension == 1)
    {
      for (int ip = 0; ip < np_local; ip++)
	{
	  L[ip].setZero();
	  for (int j = 0; j < numneigh_pn[ip]; j++)
	    {
	      in = neigh_pn[ip][j];
	      L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	    }
	}
    }
  else if ((domain->dimension == 2) && (domain->axisymmetric == false))
    {
      for (int ip = 0; ip < np_local; ip++)
	{
	  L[ip].setZero();
	  for (int j = 0; j < numneigh_pn[ip]; j++)
	    {
	      in = neigh_pn[ip][j];
	      L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	      L[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	      L[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	      L[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
	    }
	}
    }
  else if ((domain->dimension == 2) && (domain->axisymmetric == true))
    {
      for (int ip = 0; ip < np_local; ip++)
	{
	  L[ip].setZero();
	  for (int j = 0; j < numneigh_pn[ip]; j++)
	    {
	      in = neigh_pn[ip][j];
	      L[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
	      L[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
	      L[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
	      L[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
	      L[ip](2, 2) += (*vn)[in][0] * wf_pn[ip][j] / x[ip][0];
	    }
	}
    }
  else if (domain->dimension == 3)
    {
      for (int ip = 0; ip < np_local; ip++)
	{
	  L[ip].setZero();
	  for (int j = 0; j < numneigh_pn[ip]; j++)
	    {
	      in = neigh_pn[ip][j];
	      L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	      L[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	      L[ip](0,2) += (*vn)[in][0]*wfd_pn[ip][j][2];
	      L[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	      L[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
	      L[ip](1,2) += (*vn)[in][1]*wfd_pn[ip][j][2];
	      L[ip](2,0) += (*vn)[in][2]*wfd_pn[ip][j][0];
	      L[ip](2,1) += (*vn)[in][2]*wfd_pn[ip][j][1];
	      L[ip](2,2) += (*vn)[in][2]*wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_USL()
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *vn = &grid->v_update;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
      }
    }
  }
  else if (domain->dimension == 2)
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        L[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
      }
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        L[ip](0, 2) += (*vn)[in][0] * wfd_pn[ip][j][2];
        L[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
        L[ip](1, 2) += (*vn)[in][1] * wfd_pn[ip][j][2];
        L[ip](2, 0) += (*vn)[in][2] * wfd_pn[ip][j][0];
        L[ip](2, 1) += (*vn)[in][2] * wfd_pn[ip][j][1];
        L[ip](2, 2) += (*vn)[in][2] * wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_deformation_gradient()
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *xn = &grid->x;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  Eigen::Vector3d dx;
  Eigen::Matrix3d Ftemp, eye;
  eye.setIdentity();

  if (domain->dimension == 1)
    {
      for (int ip = 0; ip < np_local; ip++)
	{
	  Ftemp.setZero();
	  for (int j = 0; j < numneigh_pn[ip]; j++)
	    {
	      in = neigh_pn[ip][j];
	      dx = (*xn)[in] - (*x0n)[in];
	      Ftemp(0, 0) += dx[0] * wfd_pn[ip][j][0];
      }
      F[ip](0, 0) = Ftemp(0, 0) + 1;
    }
  }
  else if (domain->dimension == 2)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
	dx = (*xn)[in] - (*x0n)[in];
        Ftemp(0, 0) += dx[0] * wfd_pn[ip][j][0];
        Ftemp(0, 1) += dx[0] * wfd_pn[ip][j][1];
        Ftemp(1, 0) += dx[1] * wfd_pn[ip][j][0];
        Ftemp(1, 1) += dx[1] * wfd_pn[ip][j][1];
      }
      F[ip] = Ftemp + eye;
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
	dx = (*xn)[in] - (*x0n)[in];
        Ftemp(0, 0) += dx[0] * wfd_pn[ip][j][0];
        Ftemp(0, 1) += dx[0] * wfd_pn[ip][j][1];
        Ftemp(0, 2) += dx[0] * wfd_pn[ip][j][2];
        Ftemp(1, 0) += dx[1] * wfd_pn[ip][j][0];
        Ftemp(1, 1) += dx[1] * wfd_pn[ip][j][1];
        Ftemp(1, 2) += dx[1] * wfd_pn[ip][j][2];
        Ftemp(2, 0) += dx[2] * wfd_pn[ip][j][0];
        Ftemp(2, 1) += dx[2] * wfd_pn[ip][j][1];
        Ftemp(2, 2) += dx[2] * wfd_pn[ip][j][2];
      }
      // F[ip].noalias() += eye;
      F[ip] = Ftemp + eye;
    }
  }
}

void Solid::compute_rate_deformation_gradient_TL_APIC()
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  //Eigen::Vector3d *vn = grid->v;
  vector<Eigen::Vector3d> *vn = &grid->v_update;
  Eigen::Vector3d dx;

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	Fdot[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	Fdot[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
	Fdot[ip](0,1) += (*vn)[in][0]*dx[1]*wf_pn[ip][j];
	Fdot[ip](1,0) += (*vn)[in][1]*dx[0]*wf_pn[ip][j];
	Fdot[ip](1,1) += (*vn)[in][1]*dx[1]*wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	Fdot[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
	Fdot[ip](0,1) += (*vn)[in][0]*dx[1]*wf_pn[ip][j];
	Fdot[ip](0,2) += (*vn)[in][0]*dx[2]*wf_pn[ip][j];
	Fdot[ip](1,0) += (*vn)[in][1]*dx[0]*wf_pn[ip][j];
	Fdot[ip](1,1) += (*vn)[in][1]*dx[1]*wf_pn[ip][j];
	Fdot[ip](1,2) += (*vn)[in][1]*dx[2]*wf_pn[ip][j];
	Fdot[ip](2,0) += (*vn)[in][2]*dx[0]*wf_pn[ip][j];
	Fdot[ip](2,1) += (*vn)[in][2]*dx[1]*wf_pn[ip][j];
	Fdot[ip](2,2) += (*vn)[in][2]*dx[2]*wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_APIC()
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  vector<Eigen::Vector3d> *vn = &grid->v_update;
  Eigen::Vector3d dx;

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	L[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	L[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
	L[ip](0,1) += (*vn)[in][0]*dx[1]*wf_pn[ip][j];
	L[ip](1,0) += (*vn)[in][1]*dx[0]*wf_pn[ip][j];
	L[ip](1,1) += (*vn)[in][1]*dx[1]*wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*x0n)[in] - x0[ip];
	L[ip](0,0) += (*vn)[in][0]*dx[0]*wf_pn[ip][j];
	L[ip](0,1) += (*vn)[in][0]*dx[1]*wf_pn[ip][j];
	L[ip](0,2) += (*vn)[in][0]*dx[2]*wf_pn[ip][j];
	L[ip](1,0) += (*vn)[in][1]*dx[0]*wf_pn[ip][j];
	L[ip](1,1) += (*vn)[in][1]*dx[1]*wf_pn[ip][j];
	L[ip](1,2) += (*vn)[in][1]*dx[2]*wf_pn[ip][j];
	L[ip](2,0) += (*vn)[in][2]*dx[0]*wf_pn[ip][j];
	L[ip](2,1) += (*vn)[in][2]*dx[1]*wf_pn[ip][j];
	L[ip](2,2) += (*vn)[in][2]*dx[2]*wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
    }
  }
}

void Solid::update_deformation_gradient()
{
  if (mat->rigid)
    return;

  bool status, lin, nh, vol_cpdi;
  Eigen::Matrix3d eye;
  eye.setIdentity();

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  if ((method_type.compare("tlcpdi") == 0 ||
       method_type.compare("ulcpdi") == 0) &&
      (update->method->style == 1))
  {
    vol_cpdi = true;
  }
  else
    vol_cpdi = false;

  for (int ip = 0; ip < np_local; ip++)
  {

    if (is_TL)
      F[ip] += update->dt * Fdot[ip];
    else
      F[ip] = (eye + update->dt * L[ip]) * F[ip];

    Finv[ip] = F[ip].inverse();

    if (vol_cpdi)
    {
      vol[ip] = 0.5 * (xpc[nc * ip + 0][0] * xpc[nc * ip + 1][1] -
                       xpc[nc * ip + 1][0] * xpc[nc * ip + 0][1] +
                       xpc[nc * ip + 1][0] * xpc[nc * ip + 2][1] -
                       xpc[nc * ip + 2][0] * xpc[nc * ip + 1][1] +
                       xpc[nc * ip + 2][0] * xpc[nc * ip + 3][1] -
                       xpc[nc * ip + 3][0] * xpc[nc * ip + 2][1] +
                       xpc[nc * ip + 3][0] * xpc[nc * ip + 0][1] -
                       xpc[nc * ip + 0][0] * xpc[nc * ip + 3][1]);
      // rho[ip] = rho0[ip];
      J[ip] = vol[ip] / vol0[ip];
    }
    else
    {
      J[ip]   = F[ip].determinant();
      vol[ip] = J[ip] * vol0[ip];
    }


    if (J[ip] <= 0.0)
    {
      cout << "Error: J[" << ip << "]<=0.0 == " << J[ip] << endl;
      cout << "F[" << ip << "]:" << endl << F[ip] << endl;
      cout << "Fdot[" << ip << "]:" << endl << Fdot[ip] << endl;
      error->one(FLERR,"");
    }
    rho[ip] = rho0[ip] / J[ip];

    if (!nh) {
      // Only done if not Neo-Hookean:

      if (is_TL) {
        status = PolDec(F[ip], R[ip]); // polar decomposition of the deformation
                                       // gradient, F = R * U

        // In TLMPM. L is computed from Fdot:
        L[ip] = Fdot[ip] * Finv[ip];
        D[ip] = 0.5 * (R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip]);

        if (!status) {
          cout << "Polar decomposition of deformation gradient failed for "
                  "particle "
               << ip << ".\n";
          cout << "F:" << endl << F[ip] << endl;
          cout << "timestep" << endl << update->ntimestep << endl;
          error->all(FLERR, "");
        }

      } else
        D[ip] = 0.5 * (L[ip] + L[ip].transpose());
    }

    // strain_increment[ip] = update->dt * D[ip];
  }
}

void Solid::update_stress()
{
  if (mat->rigid)
    return;

  max_p_wave_speed = 0;
  double flow_stress;
  Matrix3d eye, FinvT, PK1, strain_increment;
  bool lin, nh, fluid, temp;

  if (mat->type == material->constitutive_model::LINEAR)
    lin = true;
  else
    lin = false;

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  eye.setIdentity();

  if (lin) {
    for (int ip = 0; ip < np_local; ip++) {
      strain_increment = update->dt * D[ip];
      strain_el[ip] += strain_increment;
      sigma[ip] += 2 * mat->G * strain_increment +
                   mat->lambda * strain_increment.trace() * eye;

      if (is_TL) {
	vol0PK1[ip] = vol0[ip] * J[ip] *
	  (R[ip] * sigma[ip] * R[ip].transpose()) *
	  Finv[ip].transpose();
      }
    }
  } else if (nh) {
    for (int ip = 0; ip < np_local; ip++) {
      // Neo-Hookean material:
      FinvT = Finv[ip].transpose();
      PK1 = mat->G * (F[ip] - FinvT) + mat->lambda * log(J[ip]) * FinvT;
      vol0PK1[ip] = vol0[ip] * PK1;
      sigma[ip] = 1.0 / J[ip] * (F[ip] * PK1.transpose());

      strain_el[ip] =
          0.5 * (F[ip].transpose() * F[ip] - eye); // update->dt * D[ip];
    }
  } else {

    vector<double> pH(np_local, 0);
    vector<double> plastic_strain_increment(np_local, 0);
    vector<Eigen::Matrix3d> sigma_dev;
    sigma_dev.resize(np_local);
    double tav = 0;

    for (int ip = 0; ip < np_local; ip++) {

      if (mat->temp != NULL) {
        mat->eos->compute_pressure(pH[ip], ienergy[ip], J[ip], rho[ip],
                                   damage[ip], D[ip], grid->cellsize, T[ip]);
        pH[ip] += mat->temp->compute_thermal_pressure(T[ip]);

        sigma_dev[ip] = mat->strength->update_deviatoric_stress(
            sigma[ip], D[ip], plastic_strain_increment[ip],
            eff_plastic_strain[ip], eff_plastic_strain_rate[ip], damage[ip],
            T[ip]);
      } else {
        mat->eos->compute_pressure(pH[ip], ienergy[ip], J[ip], rho[ip],
                                   damage[ip], D[ip], grid->cellsize);
        sigma_dev[ip] = mat->strength->update_deviatoric_stress(
            sigma[ip], D[ip], plastic_strain_increment[ip],
            eff_plastic_strain[ip], eff_plastic_strain_rate[ip], damage[ip]);
      }

      eff_plastic_strain[ip] += plastic_strain_increment[ip];

      // // compute a characteristic time over which to average the plastic
      // strain

      tav = 1000 * grid->cellsize / mat->signal_velocity;

      eff_plastic_strain_rate[ip] -=
          eff_plastic_strain_rate[ip] * update->dt / tav;
      eff_plastic_strain_rate[ip] += plastic_strain_increment[ip] / tav;
      eff_plastic_strain_rate[ip] = MAX(0.0, eff_plastic_strain_rate[ip]);

      if (mat->damage != NULL) {
	if (update->method->temp) {
	  mat->damage->compute_damage(damage_init[ip], damage[ip], pH[ip],
				      sigma_dev[ip], eff_plastic_strain_rate[ip],
				      plastic_strain_increment[ip], T[ip]);
	} else {
	  mat->damage->compute_damage(damage_init[ip], damage[ip], pH[ip],
				      sigma_dev[ip], eff_plastic_strain_rate[ip],
				      plastic_strain_increment[ip]);
	}
      }

      if (mat->temp != NULL) {
	flow_stress = SQRT_3_OVER_2 * sigma_dev[ip].norm();
        mat->temp->compute_heat_source(T[ip], gamma[ip], flow_stress,
                                       eff_plastic_strain_rate[ip]);
        if (is_TL)
          gamma[ip] *= vol0[ip] * mat->invcp;
        else
          gamma[ip] *= vol[ip] * mat->invcp;
      }

      if (damage[ip] == 0 || pH[ip] >= 0)
	sigma[ip] = -pH[ip] * eye + sigma_dev[ip];
      else
	sigma[ip] = -pH[ip] * (1.0 - damage[ip])* eye + sigma_dev[ip];

      if (damage[ip] > 1e-10) {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye +
            sigma_dev[ip] / (mat->G * (1 - damage[ip]));
      } else {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye +
            sigma_dev[ip] / mat->G;
      }

      if (is_TL) {
	vol0PK1[ip] = vol0[ip] * J[ip] *
	  (R[ip] * sigma[ip] * R[ip].transpose()) *
	  Finv[ip].transpose();
      }
    }
  }

  double min_h_ratio = 1.0;

  for (int ip = 0; ip < np_local; ip++) {
    if (damage[ip] >= 1.0)
      continue;

    max_p_wave_speed =
        MAX(max_p_wave_speed,
            sqrt((mat->K + FOUR_THIRD * mat->G) / rho[ip]) +
                MAX(MAX(fabs(v[ip](0)), fabs(v[ip](1))), fabs(v[ip](2))));

    if (std::isnan(max_p_wave_speed)) {
      cout << "Error: max_p_wave_speed is nan with ip=" << ip
           << ", ptag[ip]=" << ptag[ip] << ", rho0[ip]=" << rho0[ip]<< ", rho[ip]=" << rho[ip]
           << ", K=" << mat->K << ", G=" << mat->G << ", J[ip]=" << J[ip]
           << endl;
      error->one(FLERR, "");
    } else if (max_p_wave_speed < 0.0) {
      cout << "Error: max_p_wave_speed= " << max_p_wave_speed
           << " with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K
           << ", G=" << mat->G << endl;
      error->one(FLERR, "");
    }

    if (is_TL) {
      EigenSolver<Matrix3d> esF(F[ip], false);
      if (esF.info()!= Success) {
	min_h_ratio = MIN(min_h_ratio,1.0);
      } else {
	min_h_ratio = MIN(min_h_ratio,fabs(esF.eigenvalues()[0].real()));
	min_h_ratio = MIN(min_h_ratio,fabs(esF.eigenvalues()[1].real()));
	min_h_ratio = MIN(min_h_ratio,fabs(esF.eigenvalues()[2].real()));
      }

      if (min_h_ratio == 0) {
	cout << "min_h_ratio == 0 with ip=" << ip
	     << "F=\n" <<  F[ip] << endl
	     << "eigenvalues of F:" << esF.eigenvalues()[0].real() << "\t" << esF.eigenvalues()[1].real() << "\t" << esF.eigenvalues()[2].real() << endl;
	cout << "esF.info()=" << esF.info() << endl;
	error->one(FLERR, "");
      }

      // dt should also be lower than the inverse of \dot{F}e_i.
      EigenSolver<Matrix3d> esFdot(Fdot[ip], false);
      if (esFdot.info()!= Success) {
	double lambda = fabs(esFdot.eigenvalues()[0].real());
	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[1].real()));
	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[2].real()));
	dtCFL = MIN(dtCFL, 0.5/lambda);
      }
    }
  }

  dtCFL = MIN(dtCFL, grid->cellsize * min_h_ratio / max_p_wave_speed);

  if (std::isnan(dtCFL))
  {
    cout << "Error: dtCFL = " << dtCFL << "\n";
    cout << "max_p_wave_speed = " << max_p_wave_speed
         << ", grid->cellsize=" << grid->cellsize << endl;
    error->one(FLERR, "");
  }
}

void Solid::compute_inertia_tensor() {

  int in;
  Eigen::Vector3d dx;
  Eigen::Matrix3d Dtemp;

  if (domain->dimension == 2) {
    for (int ip = 0; ip < np_local; ip++) {
      Dtemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        dx = grid->x0[in] - x0[ip];
        Dtemp(0, 0) += wf_pn[ip][j] * (dx[0] * dx[0]);
        Dtemp(0, 1) += wf_pn[ip][j] * (dx[0] * dx[1]);
        Dtemp(1, 1) += wf_pn[ip][j] * (dx[1] * dx[1]);
      }
      Dtemp(1, 0) = Dtemp(0, 1);
      Dtemp(2, 2) = 1;
      // if (ip==0) cout << "1 - Dtemp[" << ip << "]=\n" << Dtemp << endl;
      Di[ip] = Dtemp.inverse();
      // if (ip==0) cout << "1 - Di[" << ip << "]=\n" << Di[ip] << endl;
    }
  } else if (domain->dimension == 3) {
    for (int ip = 0; ip < np_local; ip++) {
      Dtemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        dx = grid->x0[in] - x0[ip];
        Dtemp(0, 0) += wf_pn[ip][j] * (dx[0] * dx[0]);
        Dtemp(0, 1) += wf_pn[ip][j] * (dx[0] * dx[1]);
        Dtemp(0, 2) += wf_pn[ip][j] * (dx[0] * dx[2]);
        Dtemp(1, 1) += wf_pn[ip][j] * (dx[1] * dx[1]);
        Dtemp(1, 2) += wf_pn[ip][j] * (dx[1] * dx[2]);
        Dtemp(2, 2) += wf_pn[ip][j] * (dx[2] * dx[2]);
      }
      Dtemp(1, 0) = Dtemp(0, 1);
      Dtemp(2, 1) = Dtemp(1, 2);
      Dtemp(2, 0) = Dtemp(0, 2);
      Di[ip] = Dtemp.inverse();
      // if (ip==0) cout << "1 - Di[" << ip << "]=\n" << Di[ip] << endl;
    }
  }

  // Eigen::Matrix3d eye;
  // eye.setIdentity();

  // double cellsizeSqInv = 1.0 / (grid->cellsize * grid->cellsize);

  // for (int ip = 0; ip < np_local; ip++)
  //   {
  //     if ( form_function.compare("linear") == 0)
  // 	{
  // 	  // If the form function is linear:
  // 	  if ((Di[ip](0,0) != 16.0 / 4.0 * cellsizeSqInv ) || (Di[ip](1,1) != 16.0 / 4.0 * cellsizeSqInv) || (Di[ip](2,2) != 16.0 / 4.0 * cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di[ip] << "\n and " << 4.0 * cellsizeSqInv * eye << endl;
  // 	  // Di[ip] = 16.0 / 4.0 * cellsizeSqInv * eye;
  // 	}
  //     else if (form_function.compare("quadratic-spline") == 0)
  // 	{
  // 	  // If the form function is a quadratic spline:
  // 	  if ((Di[ip](0,0) != 4.0 * cellsizeSqInv ) || (Di[ip](1,1) != 4.0 * cellsizeSqInv) || (Di[ip](2,2) != 4.0 * cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di[ip] << "\n and " << 4.0 * cellsizeSqInv * eye << endl;
  // 	}
  //     else if (form_function.compare("cubic-spline") == 0)
  // 	{
  // 	  // If the form function is a quadratic spline:
  // 	  if ((Di[ip](0,0) != 3.0 * cellsizeSqInv ) || (Di[ip](1,1) != 3.0 * cellsizeSqInv) || (Di[ip](2,2) != 3.0 * cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di[ip] << "\n and " << 3.0 * cellsizeSqInv * eye << endl;
  // 	  // If the form function is a cubic spline:
  // 	  // Di[ip] = 3.0 * cellsizeSqInv * eye;
  // 	}
  //     else if (form_function.compare("Bernstein-quadratic") == 0)
  // 	  if ((Di[ip](0,0) != 12.0 * cellsizeSqInv ) || (Di[ip](1,1) != 12.0 * cellsizeSqInv) || (Di[ip](2,2) != 12.0 * cellsizeSqInv ))
  // 	    cout << "2 - Di[" << ip << "]=\n" << Di[ip] << "\n and " << 12.0 * cellsizeSqInv * eye << endl;
  //     //Di[ip] = 12.0 * cellsizeSqInv * eye;
  //     // cout << "2 - Di[" << ip << "]=\n" << Di[ip] << endl;
  //   }
}

void Solid::copy_particle(int i, int j) {
  ptag[j]                    = ptag[i];
  x0[j]                      = x0[i];
  x[j]                       = x[i];
  v[j]                       = v[i];
  v_update[j]                = v_update[i];
  a[j]                       = a[i];
  mbp[j]                     = mbp[i];
  f[j]                       = f[i];
  vol0[j]                    = vol0[i];
  vol[j]                     = vol[i];
  rho0[j]                    = rho0[i];
  rho[j]                     = rho[i];
  mass[j]                    = mass[i];
  eff_plastic_strain[j]      = eff_plastic_strain[i];
  eff_plastic_strain_rate[j] = eff_plastic_strain_rate[i];
  damage[j]                  = damage[i];
  damage_init[j]             = damage_init[i];
  if (update->method->temp) {
    T[j]                     = T[i];
    gamma[j]                   = gamma[i];
    q[j]                       = q[i];
  }
  ienergy[j]                 = ienergy[i];
  mask[j]                    = mask[i];
  sigma[j]                   = sigma[i];
  strain_el[j]               = strain_el[i];
  vol0PK1[j]                 = vol0PK1[i];
  L[j]                       = L[i];
  F[j]                       = F[i];
  R[j]                       = R[i];
  D[j]                       = D[i];
  Finv[j]                    = Finv[i];
  Fdot[j]                    = Fdot[i];
  J[j]                       = J[i];
  if (apic)
    Di[j]                      = Di[i];

  // if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  //   {
  //     if (update->method->style == 0)
  // 	{ // CPDI-R4
  // 	  for (int id = 0; id < domain->dimension; id++)
  // 	    {
  // 	      rp0[domain->dimension * j + id] = rp0[domain->dimension * i + id];
  // 	      rp[domain->dimension * j + id]  = rp[domain->dimension * i + id];
  // 	    }
  // 	}

  //     if (update->method->style == 1)
  // 	{ // CPDI-Q4
  // 	  for (int ic = 0; ic < nc; ic++)
  // 	    {
  // 	      xpc0[nc * j + ic] = xpc0[nc * i + ic];
  // 	      xpc[nc * j + ic]  = xpc[nc * i + ic];
  // 	    }
  // 	}
  //   }
}

void Solid::pack_particle(int i, vector<double> &buf)
{
  buf.push_back(ptag[i]);

  buf.push_back(x[i](0));
  buf.push_back(x[i](1));
  buf.push_back(x[i](2));

  buf.push_back(x0[i](0));
  buf.push_back(x0[i](1));
  buf.push_back(x0[i](2));

  buf.push_back(v[i](0));
  buf.push_back(v[i](1));
  buf.push_back(v[i](2));

  buf.push_back(v_update[i](0));
  buf.push_back(v_update[i](1));
  buf.push_back(v_update[i](2));

  buf.push_back(a[i](0));
  buf.push_back(a[i](1));
  buf.push_back(a[i](2));

  buf.push_back(mbp[i](0));
  buf.push_back(mbp[i](1));
  buf.push_back(mbp[i](2));

  buf.push_back(f[i](0));
  buf.push_back(f[i](1));
  buf.push_back(f[i](2));

  buf.push_back(vol0[i]);
  buf.push_back(vol[i]);
  
  buf.push_back(rho0[i]);
  buf.push_back(rho[i]);

  buf.push_back(mass[i]);
  buf.push_back(eff_plastic_strain[i]);
  buf.push_back(eff_plastic_strain_rate[i]);
  buf.push_back(damage[i]);
  buf.push_back(damage_init[i]);
  if (update->method->temp) {
    buf.push_back(T[i]);
    buf.push_back(gamma[i]);
    buf.push_back(q[i][0]);
    buf.push_back(q[i][1]);
    buf.push_back(q[i][2]);
  }
  buf.push_back(ienergy[i]);
  buf.push_back(mask[i]);

  buf.push_back(sigma[i](0,0));
  buf.push_back(sigma[i](1,1));
  buf.push_back(sigma[i](2,2));
  buf.push_back(sigma[i](0,1));
  buf.push_back(sigma[i](0,2));
  buf.push_back(sigma[i](1,2));

  // buf.push_back(vol0PK1[i](0,0));
  // buf.push_back(vol0PK1[i](0,1));
  // buf.push_back(vol0PK1[i](0,2));
  // buf.push_back(vol0PK1[i](1,0));
  // buf.push_back(vol0PK1[i](1,1));
  // buf.push_back(vol0PK1[i](1,2));
  // buf.push_back(vol0PK1[i](2,0));
  // buf.push_back(vol0PK1[i](2,1));
  // buf.push_back(vol0PK1[i](2,2));

  // buf.push_back(L[i](0,0));
  // buf.push_back(L[i](0,1));
  // buf.push_back(L[i](0,2));
  // buf.push_back(L[i](1,0));
  // buf.push_back(L[i](1,1));
  // buf.push_back(L[i](1,2));
  // buf.push_back(L[i](2,0));
  // buf.push_back(L[i](2,1));
  // buf.push_back(L[i](2,2));

  buf.push_back(F[i](0,0));
  buf.push_back(F[i](0,1));
  buf.push_back(F[i](0,2));
  buf.push_back(F[i](1,0));
  buf.push_back(F[i](1,1));
  buf.push_back(F[i](1,2));
  buf.push_back(F[i](2,0));
  buf.push_back(F[i](2,1));
  buf.push_back(F[i](2,2));

  buf.push_back(J[i]);
}

void Solid::unpack_particle(int &i, vector<int> list, double buf[])
{
  int m;
  for (auto j: list)
    {
      m = j;

      ptag[i] = (tagint) buf[m++];

      x[i](0) = buf[m++];
      x[i](1) = buf[m++];
      x[i](2) = buf[m++];

      x0[i](0) = buf[m++];
      x0[i](1) = buf[m++];
      x0[i](2) = buf[m++];

      v[i](0) = buf[m++];
      v[i](1) = buf[m++];
      v[i](2) = buf[m++];

      v_update[i](0) = buf[m++];
      v_update[i](1) = buf[m++];
      v_update[i](2) = buf[m++];

      a[i](0) = buf[m++];
      a[i](1) = buf[m++];
      a[i](2) = buf[m++];

      mbp[i](0) = buf[m++];
      mbp[i](1) = buf[m++];
      mbp[i](2) = buf[m++];

      f[i](0) = buf[m++];
      f[i](1) = buf[m++];
      f[i](2) = buf[m++];

      vol0[i] = buf[m++];
      vol[i] = buf[m++];

      rho0[i] = buf[m++];
      rho[i] = buf[m++];

      mass[i] = buf[m++];
      eff_plastic_strain[i] = buf[m++];
      eff_plastic_strain_rate[i] = buf[m++];
      damage[i] = buf[m++];
      damage_init[i] = buf[m++];
      if (update->method->temp) {
	T[i] = buf[m++];
	gamma[i] = buf[m++];

	q[i][0] = buf[m++];
	q[i][1] = buf[m++];
	q[i][2] = buf[m++];
      }

      ienergy[i] = buf[m++];
      mask[i] = buf[m++];

      sigma[i](0,0) = buf[m++];
      sigma[i](1,1) = buf[m++];
      sigma[i](2,2) = buf[m++];
      sigma[i](0,1) = buf[m++];
      sigma[i](0,2) = buf[m++];
      sigma[i](1,2) = buf[m++];
      sigma[i](1,0) = sigma[i](0,1);
      sigma[i](2,0) = sigma[i](0,2);
      sigma[i](2,1) = sigma[i](1,2);

      // vol0PK1[i](0,0) = buf[m++];
      // vol0PK1[i](0,1) = buf[m++];
      // vol0PK1[i](0,2) = buf[m++];
      // vol0PK1[i](1,0) = buf[m++];
      // vol0PK1[i](1,1) = buf[m++];
      // vol0PK1[i](1,2) = buf[m++];
      // vol0PK1[i](2,0) = buf[m++];
      // vol0PK1[i](2,1) = buf[m++];
      // vol0PK1[i](2,2) = buf[m++];

      // L[i](0,0) = buf[m++];
      // L[i](0,1) = buf[m++];
      // L[i](0,2) = buf[m++];
      // L[i](1,0) = buf[m++];
      // L[i](1,1) = buf[m++];
      // L[i](1,2) = buf[m++];
      // L[i](2,0) = buf[m++];
      // L[i](2,1) = buf[m++];
      // L[i](2,2) = buf[m++];

      F[i](0,0) = buf[m++];
      F[i](0,1) = buf[m++];
      F[i](0,2) = buf[m++];
      F[i](1,0) = buf[m++];
      F[i](1,1) = buf[m++];
      F[i](1,2) = buf[m++];
      F[i](2,0) = buf[m++];
      F[i](2,1) = buf[m++];
      F[i](2,2) = buf[m++];

      J[i] = buf[m++];
      i++;
    }
}

void Solid::populate(vector<string> args)
{
  cout << "Solid delimitated by region ID: " << args[2] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1)
    {
      error->all(FLERR, "Error: region ID " + args[2] + " not does not exist.\n");
    }

  if (domain->created==false)
    {
      error->all(FLERR, "The domain must be created before any solids can (create_domain(...)).");
    }

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  vector<double> limits = domain->regions[iregion]->limits();

  solidlo[0] = limits[0];
  solidhi[0] = limits[1];
  solidlo[1] = limits[2];
  solidhi[1] = limits[3];
  solidlo[2] = limits[4];
  solidhi[2] = limits[5];

  solidsublo[0] = MAX(solidlo[0], sublo[0]);
  solidsublo[1] = MAX(solidlo[1], sublo[1]);
  solidsublo[2] = MAX(solidlo[2], sublo[2]);

  solidsubhi[0] = MIN(solidhi[0], subhi[0]);
  solidsubhi[1] = MIN(solidhi[1], subhi[1]);
  solidsubhi[2] = MIN(solidhi[2], subhi[2]);

  //#ifdef DEBUG
  cout << "proc " << universe->me
       << "\tsolidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2]
       << "]\t solidsubhi=["<< solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2]
       << "]\n";

//   std::vector<double> x2plot, y2plot;
//   plt::figure_size(1200, 780);
// #endif

  // Calculate total number of particles np_local:
  int nsubx, nsuby, nsubz;
  double delta;
  double hdelta;
  double Lsubx, Lsuby, Lsubz;

  double *boundlo, *boundhi;

  delta = grid->cellsize;

  //////////////////////////////////////////////////////////
  
  // if (is_TL)
  //   {
  //   // The grid will be ajusted to the solid's domain (good for TLMPM).
  //   // and we need to create the corresponding grid:
  //     grid->init(solidlo, solidhi);

  //     boundlo = solidlo;

  //     Lsubx = solidsubhi[0] - solidsublo[0];
  //     if (domain->dimension >= 2) Lsuby = solidsubhi[1] - solidsublo[1];
  //     if (domain->dimension == 3) Lsubz = solidsubhi[2] - solidsublo[2];
  //   }
  // else
  //   {
  //     // The grid is most likely bigger than the solid's domain (good for ULMPM),
  //     // so all particles created won't lie in the region, they will need to be
  //     // checked:

  //     boundlo = domain->sublo;

  //     Lsubx = domain->subhi[0] - domain->sublo[0];
  //     if (domain->dimension >= 2)
  // 	Lsuby = domain->subhi[1] - domain->sublo[1];
  //     if (domain->dimension == 3)
  // 	Lsubz = domain->subhi[2] - domain->sublo[2];
  //   }

  // nsubx = (int) (Lsubx / delta);
  // while (nsubx * delta <= Lsubx - 0.1 * delta) nsubx++;
  // nsubx++;

  // if (domain->dimension >= 2)
  //   {
  //     nsuby = (int) (Lsuby / delta);
  //   while (nsuby * delta <= Lsuby - 0.1 * delta) nsuby++;
  //   nsuby++;
  //   }
  // else
  //   {
  //     nsuby = 1;
  //   }

  // if (domain->dimension == 3)
  //   {
  //     nsubz = (int) (Lsubz/delta);
  //     while (nsubz * delta <= Lsubz - 0.1 * delta) nsubz++;
  //     nsubz++;
  //   }
  // else
  //   {
  //     nsubz = 1;
  //   }

  
  // cout << "1--proc " << universe->me << "\t nsub=["<< nsubx << "," << nsuby << "," << nsubz << "]\n";
  //////////////////////////////////////////////////////////

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (is_TL) {
      grid->init(solidlo, solidhi);
      boundlo = solidlo;
      boundhi = solidhi;
    } else {
      boundlo = domain->boxlo;
      boundhi = domain->boxhi;
    }

    double Loffsetlo[3] = {MAX(0.0, sublo[0] - boundlo[0]),
                           MAX(0.0, sublo[1] - boundlo[1]),
                           MAX(0.0, sublo[2] - boundlo[2])};
    double Loffsethi[3] = {MAX(0.0, MIN(subhi[0], boundhi[0]) - boundlo[0]),
                           MAX(0.0, MIN(subhi[1], boundhi[1]) - boundlo[1]),
                           MAX(0.0, MIN(subhi[2], boundhi[2]) - boundlo[2])};

    int noffsetlo[3] = {(int)floor(Loffsetlo[0] / delta),
                        (int)floor(Loffsetlo[1] / delta),
                        (int)floor(Loffsetlo[2] / delta)};

    int noffsethi[3] = {(int)ceil(Loffsethi[0] / delta),
                        (int)ceil(Loffsethi[1] / delta),
                        (int)ceil(Loffsethi[2] / delta)};

    cout << "1--- proc " << universe->me << " noffsetlo=[" << noffsetlo[0]
         << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
    cout << "1--- proc " << universe->me << " noffsethi=[" << noffsethi[0]
         << "," << noffsethi[1] << "," << noffsethi[2] << "]\n";

    cout << "abs=" << abs(boundlo[0] + noffsethi[0] * delta - subhi[0])<< "]\n";
    if (universe->procneigh[0][1] >= 0 &&
        abs(boundlo[0] + noffsethi[0] * delta - subhi[0]) < 1.0e-12) {
      noffsethi[0]++;
    }
    if (domain->dimension >= 2 && universe->procneigh[1][1] >= 0 &&
        abs(boundlo[1] + noffsethi[1] * delta - subhi[1]) < 1.0e-12) {
      noffsethi[1]++;
    }
    if (domain->dimension == 3 && universe->procneigh[2][1] >= 0 &&
        abs(boundlo[2] + noffsethi[2] * delta - subhi[2]) < 1.0e-12) {
      noffsethi[2]++;
    }

    cout << "2--- proc " << universe->me << " noffsethi=[" << noffsethi[0]
         << "," << noffsethi[1] << "," << noffsethi[2] << "]\n";

    nsubx = noffsethi[0] - noffsetlo[0];
    if (domain->dimension >= 2) {
      nsuby = noffsethi[1] - noffsetlo[1];
    } else {
      nsuby = 1;
    }
    if (domain->dimension >= 3) {
      nsubz = noffsethi[2] - noffsetlo[2];
    } else {
      nsubz = 1;
    }

    if (universe->procneigh[0][1] == -1) {
      while (boundlo[0] + delta * (noffsetlo[0] + nsubx - 0.5) <
             MIN(subhi[0], boundhi[0]))
        nsubx++;
    }
    if (universe->procneigh[1][1] == -1) {
      while (boundlo[1] + delta * (noffsetlo[1] + nsuby - 0.5) <
             MIN(subhi[1], boundhi[1]))
        nsuby++;
    }
    if (universe->procneigh[2][1] == -1) {
      while (boundlo[2] + delta * (noffsetlo[2] + nsubz - 0.5) <
             MIN(subhi[2], boundhi[2]))
        nsubz++;
    }
    cout << "2--proc " << universe->me << "\t nsub=[" << nsubx << "," << nsuby
         << "," << nsubz << "]\n";
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tLsub=[" << Lsubx << "," << Lsuby << "," << Lsubz << "]\t nsub=["<< nsubx << "," << nsuby << "," << nsubz << "]\n";
// #endif

  np_local = nsubx * nsuby * nsubz;

  // Create particles:

  cout << "delta = " << delta << endl;

  int l = 0;
  double vol_;

  if (domain->dimension == 1)
    vol_ = delta;
  else if (domain->dimension == 2)
    vol_ = delta * delta;
  else
    vol_ = delta * delta * delta;

  double mass_;
  if (mat->rigid)
    mass_ = 1;
  else
    mass_ = mat->rho0 * vol_;

  int np_per_cell = (int) input->parsev(args[3]);
  double xi = 0.5;
  double lp = delta;
  int nip   = 1;
  vector<double> intpoints;

  if (np_per_cell == 1)
  {
    // One particle per cell at the center:

    xi = 0.5;
    lp *= 0.5;
    nip = 1;

    intpoints = {0, 0, 0};
  }
  else if (np_per_cell == 2)
  {
    // Quadratic elements:

    if      (domain->dimension == 1) nip = 2;
    else if (domain->dimension == 2) nip = 4;
    else                             nip = 8;

    // if (is_TL && nc == 0)
    //   xi = 0.5 / sqrt(3.0);
    // else
      xi = 0.25;

    lp *= 0.25;

    intpoints = {-xi, -xi, -xi, -xi, xi, -xi, xi, -xi, -xi, xi, xi, -xi,
                 -xi, -xi, xi,  -xi, xi, xi,  xi, -xi, xi,  xi, xi, xi};
  }
  else if (np_per_cell == 3)
  {
    // Berstein elements:

    if (nc == 0)
      xi = 0.7746 / 2;
    else
      xi = 1.0 / 3.0;

    lp *= 1.0 / 6.0;
    nip = 27;

    if (domain->dimension == 1) nip = 3;
    else if (domain->dimension == 2) nip = 9;
    else nip = 27;

    intpoints = {-xi, -xi, -xi,
		 -xi, 0, -xi,
		 -xi, xi, -xi,
		 0, -xi, -xi,
		 0, 0, -xi,
		 0, xi, -xi,
		 xi, -xi, -xi,
		 xi, 0, -xi,
		 xi, xi, -xi,
		 -xi, -xi, 0,
		 -xi, 0, 0,
		 -xi, xi, 0,
		 0, -xi, 0,
		 0, 0, 0,
		 0, xi, 0,
		 xi, -xi, 0,
		 xi, 0, 0,
		 xi, xi, 0,
		 -xi, -xi, xi,
		 -xi, 0, xi,
		 -xi, xi, xi,
		 0, -xi, xi,
		 0, 0, xi,
		 0, xi, xi,
		 xi, -xi, xi,
		 xi, 0, xi,
		 xi, xi, xi};

  } else {
    lp *= 1.0 / (2 * np_per_cell);

    if (domain->dimension == 1) {
      nip = np_per_cell;
    } else if (domain->dimension == 2) {
      nip = np_per_cell * np_per_cell;
    } else {
      nip = np_per_cell * np_per_cell * np_per_cell;
    }

    double d = 1.0 / np_per_cell;

    for (int k = 0; k < np_per_cell; k++) {
      for (int i = 0; i < np_per_cell; i++) {
	for (int j = 0; j < np_per_cell; j++) {
	  intpoints.push_back((i + 0.5) * d - 0.5);
	  intpoints.push_back((j + 0.5) * d - 0.5);
	  intpoints.push_back((k + 0.5) * d - 0.5);
	}
      }
    }
  }

  np_local *= nip;
// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tnp_local=" << np_local << endl;
// #endif
  mass_ /= (double) nip;
  vol_ /= (double) nip;

  // Allocate the space in the vectors for np particles:
  grow(np_local);

  int dim = domain->dimension;
  bool r4 = false;
  bool q4 = false;

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {
    if (update->method->style == 0)
    { // CPDI-R4
      r4 = true;
    }
    if (update->method->style == 1)
    { // CPDI-Q4
      q4 = true;
    }
  }

  for (int i = 0; i < nsubx; i++)
    {
      for (int j = 0; j < nsuby; j++)
	{
	  for (int k = 0; k < nsubz; k++)
	    {
	      for (int ip = 0; ip < nip; ip++)
		{

		  if (l >= np_local)
		    {
		      cout << "Error in Solid::populate(), exceeding the allocated number of particles.\n";
		      cout << "l = " << l << endl;
		      cout << "np_local = " << np_local << endl;
		      error->all(FLERR, "");
		    }

		  x0[l][0] = x[l][0] =
		    boundlo[0] + delta*(noffsetlo[0] + i + 0.5 + intpoints[3*ip+0]);
		  x0[l][1] = x[l][1] =
		    boundlo[1] + delta*(noffsetlo[1] + j + 0.5 + intpoints[3*ip+1]);
		  if (dim == 3)
		    x0[l][2] = x[l][2] =
		      boundlo[2] + delta*(noffsetlo[2] + k + 0.5 + intpoints[3*ip+2]);
		  else
		    x0[l][2] = x[l][2] = 0;

		  // Check if the particle is inside the region:
		  if (domain->inside_subdomain(x0[l][0], x0[l][1], x0[l][2]) && domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2]) == 1)
		    {
		      // cout << "Inside\n";

		      if (r4)
			{ // CPDI-R4
			  rp0[dim * l][0] = rp[dim * l][0] = lp;
			  rp0[dim * l][1] = rp[dim * l][1] = 0;
			  rp0[dim * l][2] = rp[dim * l][2] = 0;

			  if (dim >= 2)
			    {
			      rp0[dim * l + 1][0] = rp[dim * l + 1][0] = 0;
			      rp0[dim * l + 1][1] = rp[dim * l + 1][1] = lp;
			      rp0[dim * l + 1][2] = rp[dim * l + 1][2] = 0;

			      if (dim == 3)
				{
				  rp0[dim * l + 2][0] = rp[dim * l + 1][0] = 0;
				  rp0[dim * l + 2][1] = rp[dim * l + 1][1] = 0;
				  rp0[dim * l + 2][0] = rp[dim * l + 1][0] = lp;
				}
			    }
			}
		      if (q4)
			{ // CPDI-Q4
			  xpc0[nc * l][0] = xpc[nc * l][0] = x0[l][0] - lp;

			  xpc0[nc * l + 1][0] = xpc[nc * l + 1][0] = x0[l][0] + lp;

			  if (dim >= 2)
			    {
			      xpc0[nc * l][1] = xpc[nc * l][1] = x0[l][1] - lp;
			      xpc0[nc * l + 1][1] = xpc[nc * l + 1][1] = x0[l][1] - lp;

			      xpc0[nc * l + 2][0] = xpc[nc * l + 2][0] = x0[l][0] + lp;
			      xpc0[nc * l + 2][1] = xpc[nc * l + 2][1] = x0[l][1] + lp;

			      xpc0[nc * l + 3][0] = xpc[nc * l + 3][0] = x0[l][0] - lp;
			      xpc0[nc * l + 3][1] = xpc[nc * l + 3][1] = x0[l][1] + lp;
			    }

			  if (dim == 3)
			    {
			      xpc0[nc * l][2] = xpc[nc * l][2] = x0[l][2] - lp;
			      xpc0[nc * l + 1][2] = xpc[nc * l + 1][2] = x0[l][2] - lp;
			      xpc0[nc * l + 2][2] = xpc[nc * l + 2][2] = x0[l][2] - lp;
			      xpc0[nc * l + 3][2] = xpc[nc * l + 3][2] = x0[l][2] - lp;

			      xpc0[nc * l + 4][0] = xpc[nc * l + 4][0] = x0[l][0] - lp;
			      xpc0[nc * l + 4][1] = xpc[nc * l + 4][1] = x0[l][1] - lp;
			      xpc0[nc * l + 4][2] = xpc[nc * l + 4][2] = x0[l][2] + lp;

			      xpc0[nc * l + 5][0] = xpc[nc * l + 5][0] = x0[l][0] + lp;
			      xpc0[nc * l + 5][1] = xpc[nc * l + 5][1] = x0[l][1] - lp;
			      xpc0[nc * l + 5][2] = xpc[nc * l + 5][2] = x0[l][2] + lp;

			      xpc0[nc * l + 6][0] = xpc[nc * l + 6][0] = x0[l][0] + lp;
			      xpc0[nc * l + 6][1] = xpc[nc * l + 6][1] = x0[l][1] + lp;
			      xpc0[nc * l + 6][2] = xpc[nc * l + 6][2] = x0[l][2] + lp;

			      xpc0[nc * l + 7][0] = xpc[nc * l + 7][0] = x0[l][0] - lp;
			      xpc0[nc * l + 7][1] = xpc[nc * l + 7][1] = x0[l][1] + lp;
			      xpc0[nc * l + 7][2] = xpc[nc * l + 7][2] = x0[l][2] + lp;
			    }
			}
		      l++;
		    }
		}
	    }
	}
    }

  if (np_local > l)
    {
      grow(l);
    }
  np_local = l; // Adjust np_local to account for the particles outside the domain

  tagint ptag0 = 0;

  for (int proc=0; proc<universe->nprocs; proc++){
    int np_local_bcast;
    if (proc == universe->me) {
      // Send np_local
      np_local_bcast = np_local;
    } else {
      // Receive np_local
      np_local_bcast = 0;
    }
    MPI_Bcast(&np_local_bcast,1,MPI_INT,proc,universe->uworld);
    if (universe->me > proc) ptag0 += np_local_bcast;
  }

// #ifdef DEBUG
//   cout << "proc " << universe->me << "\tptag0 = " << ptag0 << endl;
// #endif
  np_local = l; // Adjust np to account for the particles outside the domain
  cout << "np_local=" << np_local << endl;

  for (int i = 0; i < np_local; i++)
  {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    mbp[i].setZero();
    v_update[i].setZero();
    rho0[i] = rho[i] = mat->rho0;

    if (domain->axisymmetric == true) {
      mass[i] = mass_ * x0[i][0];
      vol0[i] = vol[i] = mass[i] / rho0[i];
    } else {
      mass[i] = mass_;
      vol0[i] = vol[i] = vol_;
    }

    eff_plastic_strain[i]      = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i]                  = 0;
    damage_init[i]             = 0;
    if (update->method->temp) {
      T[i]                     = T0;
      gamma[i]                 = 0;
      q[i].setZero();
    }
    ienergy[i]                 = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
    if (apic)
      Di[i].setZero();

    J[i] = 1;
    mask[i] = 1;

    ptag[i] = ptag0 + i + 1 + domain->np_total;
  }

  if (l != np_local)
  {
    cout << "Error l=" << l << " != np_local=" << np_local << endl;
    error->one(FLERR, "");
  }

  int np_local_reduced;
  MPI_Allreduce(&np_local, &np_local_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);
  np += np_local_reduced;
  domain->np_total += np;

// #ifdef DEBUG
//   vector<double> xdomain, ydomain;
//   vector<double> xsubdomain, ysubdomain;

//   xdomain.push_back(domain->boxlo[0]);
//   ydomain.push_back(domain->boxlo[1]);

//   xdomain.push_back(domain->boxhi[0]);
//   ydomain.push_back(domain->boxlo[1]);

//   xdomain.push_back(domain->boxhi[0]);
//   ydomain.push_back(domain->boxhi[1]);

//   xdomain.push_back(domain->boxlo[0]);
//   ydomain.push_back(domain->boxhi[1]);

//   xdomain.push_back(domain->boxlo[0]);
//   ydomain.push_back(domain->boxlo[1]);


//   xsubdomain.push_back(sublo[0]);
//   ysubdomain.push_back(sublo[1]);

//   xsubdomain.push_back(subhi[0]);
//   ysubdomain.push_back(sublo[1]);

//   xsubdomain.push_back(subhi[0]);
//   ysubdomain.push_back(subhi[1]);

//   xsubdomain.push_back(sublo[0]);
//   ysubdomain.push_back(subhi[1]);

//   xsubdomain.push_back(sublo[0]);
//   ysubdomain.push_back(sublo[1]);

//   plt::plot(xdomain, ydomain, "b-");
//   plt::plot(xsubdomain, ysubdomain, "r-");
//   plt::plot(x2plot, y2plot, ".");
//   plt::axis("equal");
//   plt::save("debug-proc_" + to_string(universe->me) + ".png");
//   plt::close();
// #endif
}

void Solid::update_particle_domain()
{
  int dim = domain->dimension;

  if (update->method->style == 0)
  { // CPDI-R4
    for (int ip = 0; ip < np_local; ip++)
    {
      rp[dim * ip] = F[ip] * rp0[dim * ip];
      if (dim >= 2)
        rp[dim * ip + 1] = F[ip] * rp0[dim * ip + 1];
      if (dim == 3)
        rp[dim * ip + 2] = F[ip] * rp0[dim * ip + 2];
    }
  }
  // if (update->method->style == 1) { // CPDI-Q4
  //   for (int ip=0; ip<np; ip++) {
  //     for(int ic=0; ic<nc; ic++) {
  // 	xpc[nc*ip + ic] += update->dt*v_update[ip];
  //     }
  //   }
  // }
}

void Solid::read_file(string fileName) {}

void Solid::read_mesh(string fileName)
{

  string line;

  int id;
  int length;
  int elemType;
  vector<string> splitLine;
  array<double, 3> xn;

  int nodeCount; // Count the number of nodes
  vector<array<double, 3>> nodes;

// #ifdef DEBUG
//   std::vector<double> x2plot, y2plot;
//   std::vector<double> xcplot, ycplot;
// #endif

  ifstream file(fileName, std::ios::in);

  if (!file)
  {
    cout << "Error: unable to open mesh file " << fileName << endl;
    error->one(FLERR, "");
  }

  cout << "Reading Gmsh mesh file ...\n";

  while (getline(file, line))
  {

    if (line.compare("$MeshFormat") == 0)
    {
      // Read mesh format informations:
      double version;
      file >> version;

      if (version >= 3.0)
      {
        cout << "Gmsh mesh file version >=3.0 not supported.\n";
        error->one(FLERR, "");
      }

      getline(file, line);
      if (line.compare("$EndMeshFormat") == 0)
      {
        cout << "Reading format...done!\n";
        break;
      }
      else
        cout << "Unexpected line: " << line << ". $EndMeshFormat expected!!\n";
    }

    if (line.compare("$Nodes") == 0)
    {
      cout << "Reading nodes...\n";
      // Read mesh node informations:
      file >> nodeCount;

      nodes.resize(nodeCount);

      for (int in = 0; in < nodeCount; in++)
      {
        file >> id >> xn[0] >> xn[1] >> xn[2];
        if (abs(xn[0]) < 1.0e-12)
          xn[0] = 0;
        if (abs(xn[1]) < 1.0e-12)
          xn[1] = 0;
        if (abs(xn[2]) < 1.0e-12)
          xn[2] = 0;

        if (domain->dimension == 1 && (xn[1] != 0. || xn[2] != 0.))
        {
          cout << "Error: node " << id << " has non 1D component.\n";
          error->one(FLERR, "");
        }

        if (domain->dimension == 2 && xn[2] != 0.)
        {
          cout << "Error: node " << id << " has non zero z component.\n";
          error->one(FLERR, "");
        }

        nodes[in] = xn;

        // Adjust solid bounds:
        if (xn[0] < solidlo[0])
          solidlo[0] = xn[0];
        if (xn[1] < solidlo[1])
          solidlo[1] = xn[1];
        if (xn[2] < solidlo[2])
          solidlo[2] = xn[2];

        if (xn[0] > solidhi[0])
          solidhi[0] = xn[0];
        if (xn[1] > solidhi[1])
          solidhi[1] = xn[1];
        if (xn[2] > solidhi[2])
          solidhi[2] = xn[2];
      }

      getline(file, line);
      if (line.compare("$EndNodes") == 0)
      {
        cout << "Reading nodes...done!\n";
      }
      else
        cout << "Unexpected line: " << line << ". $EndNodes expected!!\n";
    }

    if (line.compare("$Elements") == 0)
    {
      cout << "Reading elements...\n";
      file >> np; // Number of elements
      getline(file, line);

      // Allocate the space in the vectors for np particles:
      grow(np);

      for (int ie = 0; ie < np; ie++)
      {
        getline(file, line);

        boost::split(splitLine, line, boost::is_any_of("\t "));

        length = splitLine.size();

        elemType = boost::lexical_cast<int>(splitLine[1]);
        // elemType == 1: 2-node   line element
        // elemType == 3: 4-node   quadrangle
        // elemType == 4: 4-node   tetrahedra

        if (elemType == 1)
        {
          int no1 = boost::lexical_cast<int>(splitLine[5]) - 1;
          int no2 = boost::lexical_cast<int>(splitLine[6]) - 1;

          xpc0[nc * ie][0] = xpc[nc * ie][0] = nodes[no1][0];
          xpc0[nc * ie][1] = xpc[nc * ie][1] = nodes[no1][1];
          xpc0[nc * ie][2] = xpc[nc * ie][2] = nodes[no1][2];

          xpc0[nc * ie + 1][0] = xpc[nc * ie + 1][0] = nodes[no2][0];
          xpc0[nc * ie + 1][1] = xpc[nc * ie + 1][1] = nodes[no2][1];
          xpc0[nc * ie + 1][2] = xpc[nc * ie + 1][2] = nodes[no2][2];

          x0[ie][0] = x[ie][0] = 0.5 * (nodes[no1][0] + nodes[no2][0]);
          x0[ie][1] = x[ie][1] = 0.5 * (nodes[no1][1] + nodes[no2][1]);
          x0[ie][2] = x[ie][2] = 0.5 * (nodes[no1][2] + nodes[no2][2]);
        }
        else if (elemType == 3)
        {

// #ifdef DEBUG
//           xcplot.clear();
//           xcplot.resize(5, 0);
//           ycplot.clear();
//           ycplot.resize(5, 0);
// #endif

          int no1 = boost::lexical_cast<int>(splitLine[5]) - 1;
          int no2 = boost::lexical_cast<int>(splitLine[6]) - 1;
          int no3 = boost::lexical_cast<int>(splitLine[7]) - 1;
          int no4 = boost::lexical_cast<int>(splitLine[8]) - 1;

          if (method_type.compare("tlcpdi") == 0 ||
              method_type.compare("ulcpdi") == 0)
          {

            xpc0[nc * ie][0] = xpc[nc * ie][0] = nodes[no1][0];
            xpc0[nc * ie][1] = xpc[nc * ie][1] = nodes[no1][1];
            xpc0[nc * ie][2] = xpc[nc * ie][2] = nodes[no1][2];

            xpc0[nc * ie + 1][0] = xpc[nc * ie + 1][0] = nodes[no2][0];
            xpc0[nc * ie + 1][1] = xpc[nc * ie + 1][1] = nodes[no2][1];
            xpc0[nc * ie + 1][2] = xpc[nc * ie + 1][2] = nodes[no2][2];

            xpc0[nc * ie + 2][0] = xpc[nc * ie + 2][0] = nodes[no3][0];
            xpc0[nc * ie + 2][1] = xpc[nc * ie + 2][1] = nodes[no3][1];
            xpc0[nc * ie + 2][2] = xpc[nc * ie + 2][2] = nodes[no3][2];

            xpc0[nc * ie + 3][0] = xpc[nc * ie + 3][0] = nodes[no4][0];
            xpc0[nc * ie + 3][1] = xpc[nc * ie + 3][1] = nodes[no4][1];
            xpc0[nc * ie + 3][2] = xpc[nc * ie + 3][2] = nodes[no4][2];
          }

// #ifdef DEBUG
//           xcplot[0] = xpc0[nc * ie][0];
//           ycplot[0] = xpc0[nc * ie][1];
//           xcplot[1] = xpc0[nc * ie + 1][0];
//           ycplot[1] = xpc0[nc * ie + 1][1];
//           xcplot[2] = xpc0[nc * ie + 2][0];
//           ycplot[2] = xpc0[nc * ie + 2][1];
//           xcplot[3] = xpc0[nc * ie + 3][0];
//           ycplot[3] = xpc0[nc * ie + 3][1];
//           xcplot[4] = xpc0[nc * ie][0];
//           ycplot[4] = xpc0[nc * ie][1];
//           plt::plot(xcplot, ycplot, "r-");
// #endif

          x0[ie][0] = x[ie][0] = 0.25 * (nodes[no1][0] + nodes[no2][0] +
                                         nodes[no3][0] + nodes[no4][0]);
          x0[ie][1] = x[ie][1] = 0.25 * (nodes[no1][1] + nodes[no2][1] +
                                         nodes[no3][1] + nodes[no4][1]);
          x0[ie][2] = x[ie][2] = 0.25 * (nodes[no1][2] + nodes[no2][2] +
                                         nodes[no3][2] + nodes[no4][2]);

          // vol0[ie] = vol[ie] = 0.5*(xpc[nc*ie+0][0]*xpc[nc*ie+1][1] -
          // xpc[nc*ie+1][0]*xpc[nc*ie+0][1]
          // 	    + xpc[nc*ie+1][0]*xpc[nc*ie+2][1] -
          // xpc[nc*ie+2][0]*xpc[nc*ie+1][1]
          // 	    + xpc[nc*ie+2][0]*xpc[nc*ie+3][1] -
          // xpc[nc*ie+3][0]*xpc[nc*ie+2][1]
          // 	    + xpc[nc*ie+3][0]*xpc[nc*ie+0][1] -
          // xpc[nc*ie+0][0]*xpc[nc*ie+3][1]);

          vol0[ie] = vol[ie] =
              0.5 *
              (nodes[no1][0] * nodes[no2][1] - nodes[no2][0] * nodes[no1][1] +
               nodes[no2][0] * nodes[no3][1] - nodes[no3][0] * nodes[no2][1] +
               nodes[no3][0] * nodes[no4][1] - nodes[no4][0] * nodes[no3][1] +
               nodes[no4][0] * nodes[no1][1] - nodes[no1][0] * nodes[no4][1]);

// #ifdef DEBUG
//           x2plot.push_back(x0[ie][0]);
//           y2plot.push_back(x0[ie][1]);
// #endif
        }
        else if (elemType == 4)
        {

          int no1 = boost::lexical_cast<int>(splitLine[5]) - 1;
          int no2 = boost::lexical_cast<int>(splitLine[6]) - 1;
          int no3 = boost::lexical_cast<int>(splitLine[7]) - 1;
          int no4 = boost::lexical_cast<int>(splitLine[8]) - 1;

          double x1 = nodes[no1][0];
          double y1 = nodes[no1][1];
          double z1 = nodes[no1][2];
          double x2 = nodes[no2][0];
          double y2 = nodes[no2][1];
          double z2 = nodes[no2][2];
          double x3 = nodes[no3][0];
          double y3 = nodes[no3][1];
          double z3 = nodes[no3][2];
          double x4 = nodes[no4][0];
          double y4 = nodes[no4][1];
          double z4 = nodes[no4][2];

          double x21 = x2 - x1;
          double x32 = x3 - x2;
          double x43 = x4 - x3;
          double x42 = x4 - x2;
          double y23 = y2 - y3;
          double y34 = y3 - y4;
          double y12 = y1 - y2;
          double y42 = y4 - y2;
          double z34 = z3 - z4;
          double z23 = z2 - z3;
          double z12 = z1 - z2;
          double z42 = z4 - z2;

          x0[ie][0] = x[ie][0] = 0.25 * (x1 + x2 + x3 + x4);
          x0[ie][1] = x[ie][1] = 0.25 * (y1 + y2 + y3 + y4);
          x0[ie][2] = x[ie][2] = 0.25 * (z1 + z2 + z3 + z4);

          vol0[ie] = vol[ie] = (1 / 6) * (x21 * (y23 * z34 - y34 * z23) +
                                          x32 * (y34 * z12 - y12 * z34) +
                                          x43 * (y12 * z23 - y23 * z12));
        }
        else
        {
          cout << "Element type " << elemType << " not supported!!\n";
          error->one(FLERR, "");
        }
      }

      getline(file, line);
      if (line.compare("$EndElements") == 0)
      {
        cout << "Reading elements...done!\n";
        break;
      }
      else
        cout << "Unexpected line: " << line << ". $EndElements expected!!\n";
    }
  }

  cout << "np=" << np << endl;

  for (int i = 0; i < np; i++)
  {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    mbp[i].setZero();
    v_update[i].setZero();
    rho0[i] = rho[i]           = mat->rho0;
    mass[i]                    = mat->rho0 * vol0[i];
    eff_plastic_strain[i]      = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i]                  = 0;
    damage_init[i]             = 0;
    if (update->method->temp) {
      T[i]                     = T0;
      gamma[i]                 = 0;
      q[i].setZero();
    }
    ienergy[i]                 = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
    if (apic)
      Di[i].setZero();

    J[i] = 1;

    if (x0[i][0] < solidlo[0])
      solidlo[0] = x0[i][0];
    if (x0[i][1] < solidlo[1])
      solidlo[1] = x0[i][1];
    if (x0[i][2] < solidlo[2])
      solidlo[2] = x0[i][2];

    if (x0[i][0] > solidhi[0])
      solidhi[0] = x0[i][0];
    if (x0[i][1] > solidhi[1])
      solidhi[1] = x0[i][1];
    if (x0[i][2] > solidhi[2])
      solidhi[2] = x0[i][2];
  }

  if (grid->nnodes == 0)
  {
    grid->init(solidlo, solidhi);
  }

// #ifdef DEBUG
//   plt::plot(x2plot, y2plot, ".");
// #endif
}


void Solid::compute_temperature_nodes(bool reset) {
  double Ttemp, Ttemp_update;
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++) {
    if (reset) {
      grid->T[in] = 0;
    }

    if (grid->mass[in] > 0) {
      Ttemp = 0;

      for (int j = 0; j < numneigh_np[in]; j++) {
        ip = neigh_np[in][j];
        Ttemp += wf_np[in][j] * mass[ip] * T[ip];
      }
      Ttemp /= grid->mass[in];
      grid->T[in] += Ttemp;
    }
  }
}

void Solid::compute_external_temperature_driving_forces_nodes(bool reset) {
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++) {
    if (reset)
      grid->Qext[in] = 0;

    if (grid->mass[in] > 0) {
      for (int j = 0; j < numneigh_np[in]; j++) {
        ip = neigh_np[in][j];
        grid->Qext[in] += wf_np[in][j] * gamma[ip];
      }
    }
  }
}

void Solid::compute_internal_temperature_driving_forces_nodes() {
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++) {
    grid->Qint[in] = 0;
    for (int j = 0; j < numneigh_np[in]; j++) {
      ip = neigh_np[in][j];
      grid->Qint[in] += wfd_np[in][j].dot(q[ip]);

      if (domain->axisymmetric == true) {
	error->one(FLERR,"Temperature and axisymmetric not yet supported.\n");
        //ftemp[0] -= vol0PK1[ip](2, 2) * wf_np[in][j] / x0[ip][0];
      }
    }
  }
}

void Solid::update_particle_temperature() {
  int in;
  for (int ip = 0; ip < np_local; ip++) {
    for (int j = 0; j < numneigh_pn[ip]; j++) {
      in = neigh_pn[ip][j];
      T[ip] += wf_pn[ip][j] * (grid->T_update[in] - grid->T[in]);
      // if (ptag[ip] == 3198 || ptag[ip] == 3102) {
      // 	cout << "T[" << ptag[ip] << "]=" << T[ip] << " in=" << grid->ntag[in] << " T_update[in]=" << grid->T_update[in] << " T[in]=" << grid->T[in] << endl;
      // }
    }
  }
}

void Solid::update_heat_flux() {
  int in;

  if (is_TL) {
    for (int ip = 0; ip < np_local; ip++) {
      q[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        q[ip] -= wfd_pn[ip][j] * grid->T[in];
      }
      q[ip] *= vol0[ip] * mat->invcp * mat->kappa;
    }
  } else {
    for (int ip = 0; ip < np_local; ip++) {
      q[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++) {
        in = neigh_pn[ip][j];
        q[ip] -= wfd_pn[ip][j] * grid->T[in];
      }
      q[ip] *= vol[ip] * mat->invcp * mat->kappa;
    }
  }
}
