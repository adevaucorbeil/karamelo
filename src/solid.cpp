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
#include <math.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <numeric>

#include <algorithm>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)
#define FOUR_THIRD 1.333333333333333333333333333333333333333

vector<string> split(string s, string delimiter)
{
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  string token;
  vector<string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != string::npos)
  {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }

  res.push_back(s.substr(pos_start));
  return res;
}

Solid::Solid(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == nullptr)
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

  id = args[0];

  if (universe->me == 0)
    cout << "Creating new solid with ID: " << id << endl;

  method_type = update->method_type;

  np = 0;

  if (update->method->is_CPDI)
  {
    nc = pow(2, domain->dimension);
  }
  else
    nc = 0;

  mat = nullptr;

  if (update->method->is_TL)
  {
    is_TL = true;
    grid = new Grid(mpm);
  }
  else
  {
    is_TL = false;
    grid = domain->grid;
  }

  if (update->sub_method_type == Update::SubMethodType::APIC ||
      update->sub_method_type == Update::SubMethodType::AFLIP ||
      update->sub_method_type == Update::SubMethodType::ASFLIP)
  {
    apic = true;
  }
  else
  {
    apic = false;
  }

  dtCFL = 1.0e22;
  max_p_wave_speed = 0;
  vtot = 0;
  mtot = 0;
  comm_n = 51; // Number of double to pack for particle exchange between CPUs.

  if (args[1].compare("restart") == 0)
  {
    // If the keyword restart, we are expecting to have read_restart()
    // launched right after.
    return;
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
    error->all(FLERR, "Error: not enough arguments.\n" + usage.find(args[1])->second);
  }

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

  if (update->method->temp)
    comm_n = 55; // Number of double to pack for particle exchange between CPUs.
  else
    comm_n = 50;
}

Solid::~Solid()
{
  if (is_TL)
    delete grid;
}

void Solid::init()
{
  if (universe->me == 0)
  {
    cout << "Bounds for " << id << ":\n";
    cout << "xlo xhi: " << solidlo[0] << " " << solidhi[0] << endl;
    cout << "ylo yhi: " << solidlo[1] << " " << solidhi[1] << endl;
    cout << "zlo zhi: " << solidlo[2] << " " << solidhi[2] << endl;
  }

  // Calculate total volume:
  double vtot_local = 0;
  double mtot_local = 0;
  for (int ip = 0; ip < np_local; ip++)
  {
    vtot_local += vol[ip];
    mtot_local += mass[ip];
  }

  MPI_Allreduce(&vtot_local, &vtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&mtot_local, &mtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

  if (universe->me == 0)
  {
    cout << "Solid " << id << " total volume = " << vtot << endl;
    cout << "Solid " << id << " total mass = " << mtot << endl;
  }

  if (grid->nnodes == 0)
    grid->init(solidlo, solidhi);

  if (np == 0)
  {
    error->one(FLERR, "Error: solid does not have any particles.\n");
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  // cout << "In solid::options()" << endl;
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
      error->all(FLERR, "\n");
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

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {

    if (update->method->style == 0)
    { // CPDI-R4
      rp0.resize(domain->dimension * nparticles);
      rp.resize(domain->dimension * nparticles);
    }
    if (update->method->style == 1)
    { // CPDI-Q4
      xpc0.resize(nc * nparticles);
      xpc.resize(nc * nparticles);
    }
  }

  if (method_type.compare("tlcpdi2") == 0 || method_type.compare("ulcpdi2") == 0)
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

  is_surf.resize(nparticles);

  bigint nnodes = grid->nnodes_local + grid->nnodes_ghost;

  numneigh_np.resize(nnodes);
  neigh_np.resize(nnodes);
  wf_np.resize(nnodes);
  wfd_np.resize(nnodes);

  if (mat->cp != 0)
  {
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
    if (reset)
      grid->mass[in] = 0;

    if (grid->rigid[in] && !mat->rigid)
      continue;

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
  // double mass_rigid;
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
    {
      grid->v[in].setZero();
      // grid->v_update[in].setZero();
      if (grid->rigid[in])
      {
        grid->mb[in].setZero();
      }
    }

    if (grid->rigid[in] && !mat->rigid)
      continue;

    if (grid->mass[in] > 0)
    {
      vtemp.setZero();
      if (grid->rigid[in])
        vtemp_update.setZero();

      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        if (grid->rigid[in])
        {
          vtemp_update += (wf_np[in][j] * mass[ip]) * v_update[ip];
        }
        if (update->method->ge)
        {
          vtemp += (wf_np[in][j] * mass[ip]) *
                   (v[ip] + L[ip] * (grid->x0[in] - x[ip]));
        }
        else
        {
          vtemp += wf_np[in][j] * mass[ip] * v[ip];
        }
        // grid->v[in] += (wf_np[in][j] * mass[ip]) * v[ip]/ grid->mass[in];
      }
      vtemp /= grid->mass[in];
      grid->v[in] += vtemp;
      if (grid->rigid[in])
      {
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
  Eigen::Vector3d vtemp;

  vector<Eigen::Vector3d> *pos;
  vector<Eigen::Matrix3d> *C;

  if (is_TL)
  {
    pos = &x0;
    C = &Fdot;
  }
  else
  {
    pos = &x;
    C = &L;
  }

  for (int in = 0; in < nn; in++)
  {
    if (reset)
      grid->v[in].setZero();

    if (grid->rigid[in] && !mat->rigid)
      continue;

    if (grid->mass[in] > 0)
    {
      vtemp.setZero();
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        vtemp += (wf_np[in][j] * mass[ip]) *
                 (v[ip] + (*C)[ip] * (grid->x0[in] - (*pos)[ip]));
      }
      vtemp /= grid->mass[in];
      grid->v[in] += vtemp;
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

void Solid::compute_external_and_internal_forces_nodes_UL(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
    {
      grid->f[in].setZero();
      grid->mb[in].setZero();
    }

    if (grid->rigid[in])
    {
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
          grid->f[in][0] -=
              vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
        }
      }
    }
    else
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        grid->f[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
        grid->mb[in] += wf_np[in][j] * mbp[ip];
      }

      if (domain->axisymmetric == true)
      {
        for (int j = 0; j < numneigh_np[in]; j++)
        {
          ip = neigh_np[in][j];
          grid->f[in][0] -=
              vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
        }
      }
    }
  }
}

void Solid::compute_external_and_internal_forces_nodes_UL_MLS(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;
  vector<Eigen::Vector3d> *pos;

  if (is_TL)
  {
    pos = &x0;
  }
  else
  {
    pos = &x;
  }

  for (int in = 0; in < nn; in++)
  {
    if (reset)
    {
      grid->f[in].setZero();
      grid->mb[in].setZero();
    }

    if (grid->rigid[in])
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        // grid->f[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
        grid->f[in] -= vol[ip] * wf_np[in][j] *
                       (sigma[ip] * Di * (grid->x0[in] - (*pos)[ip]));
      }

      if (domain->axisymmetric == true)
      {
        for (int j = 0; j < numneigh_np[in]; j++)
        {
          ip = neigh_np[in][j];
          grid->f[in][0] -=
              vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
        }
      }
    }
    else
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        // grid->f[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
        grid->f[in] -= vol[ip] * wf_np[in][j] *
                       (sigma[ip] * Di * (grid->x0[in] - (*pos)[ip]));
        grid->mb[in] += wf_np[in][j] * mbp[ip];
      }

      if (domain->axisymmetric == true)
      {
        for (int j = 0; j < numneigh_np[in]; j++)
        {
          ip = neigh_np[in][j];
          grid->f[in][0] -=
              vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
        }
      }
    }
  }
}

void Solid::compute_particle_accelerations_velocities_and_positions()
{

  vector<Eigen::Vector3d> vc_update;
  vc_update.resize(nc);

  int in;

  bool update_corners;
  double inv_dt = 1.0 / update->dt;
  Vector3d dummy;

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
    a[ip].setZero();
    if (update_corners)
      for (int i = 0; i < nc; i++)
        vc_update[i].setZero();

    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * grid->v_update[in];
      a[ip] += wf_pn[ip][j] * (grid->v_update[in] - grid->v[in]);
      // x[ip] += update->dt * wf_pn[ip][j] * grid->v_update[in];

      if (update_corners)
      {
        for (int ic = 0; ic < nc; ic++)
        {
          vc_update[ic] += wf_pn_corners[nc * ip + ic][j] * grid->v_update[in];
        }
      }
    }
    a[ip] *= inv_dt;
    f[ip] = a[ip] * mass[ip];
    x[ip] += update->dt * v_update[ip];

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

void Solid::compute_particle_accelerations_velocities()
{

  vector<Eigen::Vector3d> vc_update;
  vc_update.resize(nc);

  int in;

  bool update_corners;
  double inv_dt = 1.0 / update->dt;
  Vector3d dummy;

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
    a[ip].setZero();
    if (update_corners)
      for (int i = 0; i < nc; i++)
        vc_update[i].setZero();

    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * grid->v_update[in];
      a[ip] += wf_pn[ip][j] * (grid->v_update[in] - grid->v[in]);

      if (update_corners)
      {
        for (int ic = 0; ic < nc; ic++)
        {
          vc_update[ic] += wf_pn_corners[nc * ip + ic][j] * grid->v_update[in];
        }
      }
    }
    a[ip] *= inv_dt;
    f[ip] = a[ip] * mass[ip];

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
  double inv_dt = 1.0 / update->dt;

  int in;

  for (int ip = 0; ip < np_local; ip++)
  {
    a[ip].setZero();
    if (mat->rigid)
      continue;
    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      a[ip] += wf_pn[ip][j] * (grid->v_update[in] - grid->v[in]);
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
  }
}
void Solid::update_particle_velocities_and_positions(double alpha)
{
  for (int ip = 0; ip < np_local; ip++)
  {
    v[ip] = (1 - alpha) * v_update[ip] + alpha * (v[ip] + update->dt * a[ip]);
    x[ip] += update->dt * v[ip];
  }
}

void Solid::compute_rate_deformation_gradient_TL(bool doublemapping)
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *vn;

  if (doublemapping)
    vn = &grid->v;
  else
    vn = &grid->v_update;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
      }
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == true))
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
        Fdot[ip](2, 2) += (*vn)[in][0] * wf_pn[ip][j] / x0[ip][0];
      }
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == false))
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += (*vn)[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += (*vn)[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](1, 0) += (*vn)[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += (*vn)[in][1] * wfd_pn[ip][j][1];
      }
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
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

void Solid::compute_rate_deformation_gradient_UL(bool doublemapping)
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *vn;

  if (doublemapping)
    vn = &grid->v;
  else
    vn = &grid->v_update;

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
  else if ((domain->dimension == 2) && (domain->axisymmetric == false))
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

void Solid::compute_rate_deformation_gradient_TL_APIC(bool doublemapping)
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *x0n = &grid->x0;

  vector<Eigen::Vector3d> *vn;

  if (doublemapping)
    vn = &grid->v;
  else
    vn = &grid->v_update;

  Eigen::Vector3d dx;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = (*x0n)[in] - x0[ip];
        Fdot[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
      }
      Fdot[ip] *= Di;
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == false))
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = (*x0n)[in] - x0[ip];
        Fdot[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        Fdot[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        Fdot[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        Fdot[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
      }
      Fdot[ip] *= Di;
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == true))
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = (*x0n)[in] - x0[ip];
        Fdot[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        Fdot[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        Fdot[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        Fdot[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
        Fdot[ip](2, 2) += (*vn)[in][0] * wf_pn[ip][j] / x0[ip][0];
      }
      Fdot[ip] *= Di;
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = (*x0n)[in] - x0[ip];
        Fdot[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        Fdot[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        Fdot[ip](0, 2) += (*vn)[in][0] * dx[2] * wf_pn[ip][j];
        Fdot[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        Fdot[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
        Fdot[ip](1, 2) += (*vn)[in][1] * dx[2] * wf_pn[ip][j];
        Fdot[ip](2, 0) += (*vn)[in][2] * dx[0] * wf_pn[ip][j];
        Fdot[ip](2, 1) += (*vn)[in][2] * dx[1] * wf_pn[ip][j];
        Fdot[ip](2, 2) += (*vn)[in][2] * dx[2] * wf_pn[ip][j];
      }
      Fdot[ip] *= Di;
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_APIC(bool doublemapping)
{
  if (mat->rigid)
    return;

  int in;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  vector<Eigen::Vector3d> *vn;

  if (doublemapping)
    vn = &grid->v;
  else
    vn = &grid->v_update;

  Eigen::Vector3d dx;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = (*x0n)[in] - x[ip];
        L[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
      }
      L[ip] *= Di;
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
        dx = (*x0n)[in] - x[ip];
        L[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        L[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        L[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        L[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
      }
      L[ip] *= Di;
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
        dx = (*x0n)[in] - x[ip];
        L[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        L[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        L[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        L[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
        L[ip](2, 2) += (*vn)[in][0] * wf_pn[ip][j] / x[ip][0];
      }
      L[ip] *= Di;
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
        dx = (*x0n)[in] - x[ip];
        L[ip](0, 0) += (*vn)[in][0] * dx[0] * wf_pn[ip][j];
        L[ip](0, 1) += (*vn)[in][0] * dx[1] * wf_pn[ip][j];
        L[ip](0, 2) += (*vn)[in][0] * dx[2] * wf_pn[ip][j];
        L[ip](1, 0) += (*vn)[in][1] * dx[0] * wf_pn[ip][j];
        L[ip](1, 1) += (*vn)[in][1] * dx[1] * wf_pn[ip][j];
        L[ip](1, 2) += (*vn)[in][1] * dx[2] * wf_pn[ip][j];
        L[ip](2, 0) += (*vn)[in][2] * dx[0] * wf_pn[ip][j];
        L[ip](2, 1) += (*vn)[in][2] * dx[1] * wf_pn[ip][j];
        L[ip](2, 2) += (*vn)[in][2] * dx[2] * wf_pn[ip][j];
      }
      L[ip] *= Di;
    }
  }
}

void Solid::update_deformation_gradient()
{
  if (mat->rigid)
    return;

  bool status, nh, vol_cpdi;
  Eigen::Matrix3d eye;
  eye.setIdentity();

  nh = mat->type == material->constitutive_model::NEO_HOOKEAN;

  vol_cpdi = (!method_type.compare("tlcpdi") ||
              !method_type.compare("ulcpdi")) &&
             (update->method->style == 1);

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
      J[ip] = F[ip].determinant();
      vol[ip] = J[ip] * vol0[ip];
    }

    if (J[ip] <= 0.0 && damage[ip] < 1.0)
    {
      cout << "Error: J[" << ptag[ip] << "]<=0.0 == " << J[ip] << endl;
      cout << "F[" << ptag[ip] << "]:" << endl
           << F[ip] << endl;
      cout << "Fdot[" << ptag[ip] << "]:" << endl
           << Fdot[ip] << endl;
      cout << "damage[" << ptag[ip] << "]:" << endl
           << damage[ip] << endl;
      cout << "np local: " << np_local << ", np total: " << np << " " << endl;
      error->one(FLERR, "");
    }
    rho[ip] = rho0[ip] / J[ip];

    if (!nh)
    {
      // Only done if not Neo-Hookean:

      if (is_TL)
      {
        status = PolDec(F[ip], R[ip]); // polar decomposition of the deformation
                                       // gradient, F = R * U

        // In TLMPM. L is computed from Fdot:
        L[ip] = Fdot[ip] * Finv[ip];
        D[ip] = 0.5 * (R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip]);

        if (!status)
        {
          cout << "Polar decomposition of deformation gradient failed for "
                  "particle "
               << ip << ".\n";
          cout << "F:" << endl
               << F[ip] << endl;
          cout << "timestep" << endl
               << update->ntimestep << endl;
          error->one(FLERR, "");
        }
      }
      else
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
  bool lin, nh;

  if (mat->type == material->constitutive_model::LINEAR)
    lin = true;
  else
    lin = false;

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  eye.setIdentity();

  if (lin)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      strain_increment = update->dt * D[ip];
      strain_el[ip] += strain_increment;
      sigma[ip] += 2 * mat->G * strain_increment +
                   mat->lambda * strain_increment.trace() * eye;

      if (is_TL)
      {
        vol0PK1[ip] = vol0[ip] * J[ip] *
                      (R[ip] * sigma[ip] * R[ip].transpose()) *
                      Finv[ip].transpose();
      }
    }
  }
  else if (nh)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      // Neo-Hookean material:
      FinvT = Finv[ip].transpose();
      PK1 = mat->G * (F[ip] - FinvT) + mat->lambda * log(J[ip]) * FinvT;
      vol0PK1[ip] = vol0[ip] * PK1;
      sigma[ip] = 1.0 / J[ip] * (F[ip] * PK1.transpose());

      strain_el[ip] =
          0.5 * (F[ip].transpose() * F[ip] - eye); // update->dt * D[ip];
    }
  }
  else
  {

    vector<double> pH(np_local, 0);
    vector<double> plastic_strain_increment(np_local, 0);
    vector<Eigen::Matrix3d> sigma_dev;
    sigma_dev.resize(np_local);
    double tav = 0;

    for (int ip = 0; ip < np_local; ip++)
    {

      if (mat->cp != 0)
      {
        mat->eos->compute_pressure(pH[ip], ienergy[ip], J[ip], rho[ip],
                                   damage[ip], D[ip], grid->cellsize, T[ip]);
        pH[ip] += mat->temp->compute_thermal_pressure(T[ip]);

        sigma_dev[ip] = mat->strength->update_deviatoric_stress(
            sigma[ip], D[ip], plastic_strain_increment[ip],
            eff_plastic_strain[ip], eff_plastic_strain_rate[ip], damage[ip],
            T[ip]);
      }
      else
      {
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

      if (mat->damage != nullptr)
      {
        if (update->method->temp)
        {
          mat->damage->compute_damage(damage_init[ip], damage[ip], pH[ip],
                                      sigma_dev[ip], eff_plastic_strain_rate[ip],
                                      plastic_strain_increment[ip], T[ip]);
        }
        else
        {
          mat->damage->compute_damage(damage_init[ip], damage[ip], pH[ip],
                                      sigma_dev[ip], eff_plastic_strain_rate[ip],
                                      plastic_strain_increment[ip]);
        }
      }

      if (mat->cp != 0)
      {
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
        sigma[ip] = -pH[ip] * (1.0 - damage[ip]) * eye + sigma_dev[ip];

      if (damage[ip] > 1e-10)
      {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye +
            sigma_dev[ip] / (mat->G * (1 - damage[ip]));
      }
      else
      {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye +
            sigma_dev[ip] / mat->G;
      }

      if (is_TL)
      {
        vol0PK1[ip] = vol0[ip] * J[ip] *
                      (R[ip] * sigma[ip] * R[ip].transpose()) *
                      Finv[ip].transpose();
      }
    }
  }

  double min_h_ratio = 1.0;

  for (int ip = 0; ip < np_local; ip++)
  {
    if (damage[ip] >= 1.0)
      continue;

    max_p_wave_speed =
        MAX(max_p_wave_speed,
            sqrt((mat->K + FOUR_THIRD * mat->G) / rho[ip]) +
                MAX(MAX(fabs(v[ip](0)), fabs(v[ip](1))), fabs(v[ip](2))));

    if (std::isnan(max_p_wave_speed))
    {
      cout << "Error: max_p_wave_speed is nan with ip=" << ip
           << ", ptag[ip]=" << ptag[ip] << ", rho0[ip]=" << rho0[ip] << ", rho[ip]=" << rho[ip]
           << ", K=" << mat->K << ", G=" << mat->G << ", J[ip]=" << J[ip]
           << endl;
      error->one(FLERR, "");
    }
    else if (max_p_wave_speed < 0.0)
    {
      cout << "Error: max_p_wave_speed= " << max_p_wave_speed
           << " with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K
           << ", G=" << mat->G << endl;
      error->one(FLERR, "");
    }

    if (is_TL)
    {
      EigenSolver<Matrix3d> esF(F[ip], false);
      if (esF.info() != Success)
      {
        min_h_ratio = MIN(min_h_ratio, 1.0);
      }
      else
      {
        min_h_ratio = MIN(min_h_ratio, fabs(esF.eigenvalues()[0].real()));
        min_h_ratio = MIN(min_h_ratio, fabs(esF.eigenvalues()[1].real()));
        min_h_ratio = MIN(min_h_ratio, fabs(esF.eigenvalues()[2].real()));
      }

      if (min_h_ratio == 0)
      {
        cout << "min_h_ratio == 0 with ip=" << ip
             << "F=\n"
             << F[ip] << endl
             << "eigenvalues of F:" << esF.eigenvalues()[0].real() << "\t" << esF.eigenvalues()[1].real() << "\t" << esF.eigenvalues()[2].real() << endl;
        cout << "esF.info()=" << esF.info() << endl;
        error->one(FLERR, "");
      }

      // // dt should also be lower than the inverse of \dot{F}e_i.
      // EigenSolver<Matrix3d> esFdot(Fdot[ip], false);
      // if (esFdot.info()!= Success) {
      // 	double lambda = fabs(esFdot.eigenvalues()[0].real());
      // 	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[1].real()));
      // 	lambda = MAX(lambda, fabs(esFdot.eigenvalues()[2].real()));
      // 	dtCFL = MIN(dtCFL, 0.5/lambda);
      // }
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

void Solid::compute_inertia_tensor()
{
  Eigen::Vector3d dx;

  vector<Eigen::Vector3d> *pos;
  Eigen::Matrix3d eye, Dtemp;
  eye.setIdentity();
  double cellsizeSqInv = 1.0 / (grid->cellsize * grid->cellsize);

  if (update->shape_function == Update::ShapeFunctions::LINEAR)
  {
    if (is_TL)
      pos = &x0;
    else
    {
      error->all(FLERR, "Shape function not supported for APIC and ULMPM.\n");
    }
    if (np_per_cell == 1)
    {
      Di = 16.0 / 4.0 * cellsizeSqInv * eye;
    }
    else if (np_per_cell == 2)
    {
      Di = 16.0 / 3.0 * cellsizeSqInv * eye;
    }
    else
    {
      error->all(FLERR, "Number of particle per cell not supported with linear "
                        "shape functions and APIC.\n");
    }
  }
  else if (update->shape_function == Update::ShapeFunctions::CUBIC_SPLINE)
  {

    Di = 3.0 * cellsizeSqInv * eye;
  }
  else if (update->shape_function ==
           Update::ShapeFunctions::QUADRATIC_SPLINE)
  {

    Di = 4.0 * cellsizeSqInv * eye;
  }
  else
  {
    error->all(FLERR, "Shape function not supported for APIC.\n");
  }

  if (domain->dimension == 1)
  {
    Di(1, 1) = 1;
    Di(2, 2) = 1;
  }
  else if (domain->dimension == 2)
  {
    Di(2, 2) = 1;
  }

  // if (is_TL)
  //   pos = &x0;
  // else
  //   pos = &x;

  // if (domain->dimension == 2) {
  //   for (int ip = 0; ip < np_local; ip++) {
  //     Dtemp.setZero();
  //     for (int j = 0; j < numneigh_pn[ip]; j++) {
  //       in = neigh_pn[ip][j];
  //       dx = grid->x0[in] - (*pos)[ip];
  //       Dtemp(0, 0) += wf_pn[ip][j] * (dx[0] * dx[0]);
  //       Dtemp(0, 1) += wf_pn[ip][j] * (dx[0] * dx[1]);
  //       Dtemp(1, 1) += wf_pn[ip][j] * (dx[1] * dx[1]);
  //     }
  //     Dtemp(1, 0) = Dtemp(0, 1);
  //     Dtemp(2, 2) = 1;
  //     Di[ip] = Dtemp.inverse();
  //   }
  // } else if (domain->dimension == 3) {
  //   for (int ip = 0; ip < np_local; ip++) {
  //     Dtemp.setZero();
  //     for (int j = 0; j < numneigh_pn[ip]; j++) {
  //       in = neigh_pn[ip][j];
  //       dx = grid->x0[in] - (*pos)[ip];
  //       Dtemp(0, 0) += wf_pn[ip][j] * (dx[0] * dx[0]);
  //       Dtemp(0, 1) += wf_pn[ip][j] * (dx[0] * dx[1]);
  //       Dtemp(0, 2) += wf_pn[ip][j] * (dx[0] * dx[2]);
  //       Dtemp(1, 1) += wf_pn[ip][j] * (dx[1] * dx[1]);
  //       Dtemp(1, 2) += wf_pn[ip][j] * (dx[1] * dx[2]);
  //       Dtemp(2, 2) += wf_pn[ip][j] * (dx[2] * dx[2]);
  //     }
  //     Dtemp(1, 0) = Dtemp(0, 1);
  //     Dtemp(2, 1) = Dtemp(1, 2);
  //     Dtemp(2, 0) = Dtemp(0, 2);
  //     Di[ip] = Dtemp.inverse();
  //   }
  // }

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

void Solid::copy_particle(int i, int j)
{
  ptag[j] = ptag[i];
  x0[j] = x0[i];
  x[j] = x[i];
  v[j] = v[i];
  v_update[j] = v_update[i];
  a[j] = a[i];
  mbp[j] = mbp[i];
  f[j] = f[i];
  vol0[j] = vol0[i];
  vol[j] = vol[i];
  rho0[j] = rho0[i];
  rho[j] = rho[i];
  mass[j] = mass[i];
  eff_plastic_strain[j] = eff_plastic_strain[i];
  eff_plastic_strain_rate[j] = eff_plastic_strain_rate[i];
  damage[j] = damage[i];
  damage_init[j] = damage_init[i];
  if (update->method->temp)
  {
    T[j] = T[i];
    gamma[j] = gamma[i];
    q[j] = q[i];
  }
  ienergy[j] = ienergy[i];
  mask[j] = mask[i];
  sigma[j] = sigma[i];
  strain_el[j] = strain_el[i];
  vol0PK1[j] = vol0PK1[i];
  L[j] = L[i];
  F[j] = F[i];
  R[j] = R[i];
  D[j] = D[i];
  Finv[j] = Finv[i];
  Fdot[j] = Fdot[i];
  J[j] = J[i];
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
  if (update->method->temp)
  {
    buf.push_back(T[i]);
    buf.push_back(gamma[i]);
    buf.push_back(q[i][0]);
    buf.push_back(q[i][1]);
    buf.push_back(q[i][2]);
  }
  buf.push_back(ienergy[i]);
  buf.push_back(mask[i]);

  buf.push_back(sigma[i](0, 0));
  buf.push_back(sigma[i](1, 1));
  buf.push_back(sigma[i](2, 2));
  buf.push_back(sigma[i](0, 1));
  buf.push_back(sigma[i](0, 2));
  buf.push_back(sigma[i](1, 2));

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

  buf.push_back(F[i](0, 0));
  buf.push_back(F[i](0, 1));
  buf.push_back(F[i](0, 2));
  buf.push_back(F[i](1, 0));
  buf.push_back(F[i](1, 1));
  buf.push_back(F[i](1, 2));
  buf.push_back(F[i](2, 0));
  buf.push_back(F[i](2, 1));
  buf.push_back(F[i](2, 2));

  buf.push_back(J[i]);

  buf.push_back(is_surf[i]);
}

inline void Solid::unpack_particle(const int i, const int offset, const vector<double> &buf)
{
  int m = offset;

  ptag[i] = (tagint)buf[m++];

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
  if (update->method->temp)
  {
    T[i] = buf[m++];
    gamma[i] = buf[m++];

    q[i][0] = buf[m++];
    q[i][1] = buf[m++];
    q[i][2] = buf[m++];
  }

  ienergy[i] = buf[m++];
  mask[i] = buf[m++];

  sigma[i](0, 0) = buf[m++];
  sigma[i](1, 1) = buf[m++];
  sigma[i](2, 2) = buf[m++];
  sigma[i](0, 1) = buf[m++];
  sigma[i](0, 2) = buf[m++];
  sigma[i](1, 2) = buf[m++];
  sigma[i](1, 0) = sigma[i](0, 1);
  sigma[i](2, 0) = sigma[i](0, 2);
  sigma[i](2, 1) = sigma[i](1, 2);

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

  F[i](0, 0) = buf[m++];
  F[i](0, 1) = buf[m++];
  F[i](0, 2) = buf[m++];
  F[i](1, 0) = buf[m++];
  F[i](1, 1) = buf[m++];
  F[i](1, 2) = buf[m++];
  F[i](2, 0) = buf[m++];
  F[i](2, 1) = buf[m++];
  F[i](2, 2) = buf[m++];

  J[i] = buf[m++];

  is_surf[i] = (int)buf[m++];
}

void Solid::unpack_particle(int &i, vector<int> list, vector<double> &buf)
{
  int m;
  for (auto j : list)
  {
    m = j;
    unpack_particle(i, j, buf);
    i++;
  }
}

void Solid::populate(vector<string> args)
{
  if (universe->me == 0)
  {
    cout << "Solid delimitated by region ID: " << args[2] << endl;
  }

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1)
  {
    error->all(FLERR, "Error: region ID " + args[2] + " not does not exist.\n");
  }

  if (domain->created == false)
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

#ifdef DEBUG
  cout << "proc " << universe->me
       << "\tsolidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2]
       << "]\t solidsubhi=[" << solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2]
       << "]\n";

//   std::vector<double> x2plot, y2plot;
#endif

  // Calculate total number of particles np_local:
  int nsubx, nsuby, nsubz;
  double delta;
  // double hdelta;
  // double Lsubx, Lsuby, Lsubz;

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
  if (is_TL)
  {
    grid->init(solidlo, solidhi);
    boundlo = solidlo;
    boundhi = solidhi;
  }
  else
  {
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

  // cout << "abs=" << abs(boundlo[0] + noffsethi[0] * delta - subhi[0])<< "]\n";
  if (universe->procneigh[0][1] >= 0 &&
      abs(boundlo[0] + noffsethi[0] * delta - subhi[0]) < 1.0e-12)
  {
    noffsethi[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][1] >= 0 &&
      abs(boundlo[1] + noffsethi[1] * delta - subhi[1]) < 1.0e-12)
  {
    noffsethi[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][1] >= 0 &&
      abs(boundlo[2] + noffsethi[2] * delta - subhi[2]) < 1.0e-12)
  {
    noffsethi[2]++;
  }

  cout << "2--- proc " << universe->me << " noffsethi=[" << noffsethi[0]
       << "," << noffsethi[1] << "," << noffsethi[2] << "]\n";

  nsubx = MAX(0, noffsethi[0] - noffsetlo[0]);
  if (domain->dimension >= 2)
  {
    nsuby = MAX(0, noffsethi[1] - noffsetlo[1]);
  }
  else
  {
    nsuby = 1;
  }
  if (domain->dimension >= 3)
  {
    nsubz = MAX(0, noffsethi[2] - noffsetlo[2]);
  }
  else
  {
    nsubz = 1;
  }

  if (universe->procneigh[0][1] == -1)
  {
    while (boundlo[0] + delta * (noffsetlo[0] + nsubx - 0.5) <
           MIN(subhi[0], boundhi[0]))
      nsubx++;
  }
  if (universe->procneigh[1][1] == -1)
  {
    while (boundlo[1] + delta * (noffsetlo[1] + nsuby - 0.5) <
           MIN(subhi[1], boundhi[1]))
      nsuby++;
  }
  if (universe->procneigh[2][1] == -1)
  {
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

  if (universe->me == 0)
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

  np_per_cell = (int)input->parsev(args[3]);
  double xi = 0.5;
  double lp = delta;
  int nip = 1;
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

    if (domain->dimension == 1)
      nip = 2;
    else if (domain->dimension == 2)
      nip = 4;
    else
      nip = 8;

    // if (is_TL && nc == 0)
    //   xi = 0.5 / sqrt(3.0);
    // else
    xi = 0.25;

    lp *= 0.25;

    intpoints = {-xi, -xi, -xi, -xi, xi, -xi, xi, -xi, -xi, xi, xi, -xi,
                 -xi, -xi, xi, -xi, xi, xi, xi, -xi, xi, xi, xi, xi};
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

    if (domain->dimension == 1)
      nip = 3;
    else if (domain->dimension == 2)
      nip = 9;
    else
      nip = 27;

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
  }
  else
  {
    lp *= 1.0 / (2 * np_per_cell);

    if (domain->dimension == 1)
    {
      nip = np_per_cell;
    }
    else if (domain->dimension == 2)
    {
      nip = np_per_cell * np_per_cell;
    }
    else
    {
      nip = np_per_cell * np_per_cell * np_per_cell;
    }

    double d = 1.0 / np_per_cell;

    for (int k = 0; k < np_per_cell; k++)
    {
      for (int i = 0; i < np_per_cell; i++)
      {
        for (int j = 0; j < np_per_cell; j++)
        {
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
  mass_ /= (double)nip;
  vol_ /= (double)nip;

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
              boundlo[0] + delta * (noffsetlo[0] + i + 0.5 + intpoints[3 * ip + 0]);
          x0[l][1] = x[l][1] =
              boundlo[1] + delta * (noffsetlo[1] + j + 0.5 + intpoints[3 * ip + 1]);
          if (dim == 3)
            x0[l][2] = x[l][2] =
                boundlo[2] + delta * (noffsetlo[2] + k + 0.5 + intpoints[3 * ip + 2]);
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

  for (int proc = 0; proc < universe->nprocs; proc++)
  {
    int np_local_bcast;
    if (proc == universe->me)
    {
      // Send np_local
      np_local_bcast = np_local;
    }
    else
    {
      // Receive np_local
      np_local_bcast = 0;
    }
    MPI_Bcast(&np_local_bcast, 1, MPI_INT, proc, universe->uworld);
    if (universe->me > proc)
      ptag0 += np_local_bcast;
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

    if (domain->axisymmetric == true)
    {
      mass[i] = mass_ * x0[i][0];
      vol0[i] = vol[i] = mass[i] / rho0[i];
    }
    else
    {
      mass[i] = mass_;
      vol0[i] = vol[i] = vol_;
    }

    eff_plastic_strain[i] = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
    if (update->method->temp)
    {
      T[i] = T0;
      gamma[i] = 0;
      q[i].setZero();
    }
    ienergy[i] = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
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
  domain->np_local += np_local;
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

  if (universe->me == 0)
    cout << "Reading Gmsh mesh file ...\n";

  int l = 0;

  while (getline(file, line))
  {

    if (line.compare("$MeshFormat") == 0)
    {
      // Read mesh format informations:
      double version;
      file >> version;
      cout << "version=" << version << endl;

      getline(file, line);

      // int file_type;
      // file >> file_type;
      // cout << "file_type=" << file_type << endl;

      // int data_size;
      // file >> data_size;
      // cout << "data_size=" << data_size << endl;

      if (version >= 3.0)
      {
        cout << "Gmsh mesh file version >=3.0 not supported.\n";
        error->one(FLERR, "");
      }

      getline(file, line);
      if (line.compare("$EndMeshFormat") == 0) {
	if (universe->me == 0) {
	  cout << "Reading format...done!\n";
	}
      } else {
	error->all(FLERR, "Unexpected line: " + line + ". $EndMeshFormat expected!!\n");
      }
    }

    if (line.compare("$Nodes") == 0)
    {
      if (universe->me == 0)
        cout << "Reading nodes...\n";
      // Read mesh node informations:
      file >> nodeCount;

      cout << "nodeCount=" << nodeCount << endl;

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
      getline(file, line);
      if (line.compare("$EndNodes") == 0)
      {
        if (universe->me == 0)
          cout << "Reading nodes...done!\n";
      }
      else
	error->all(FLERR, "Unexpected line: " + line + ". $EndNodes expected!!\n");
    }

    if (line.compare("$Elements") == 0)
    {
      if (universe->me == 0)
	cout << "Reading elements...\n";
      file >> np_local; // Number of elements
      getline(file, line);

      // Allocate the space in the vectors for np particles:
      grow(np_local);

      for (int ie = 0; ie < np_local; ie++)
      {
        getline(file, line);

        splitLine = split(line, " ");
	// Remove all empty elements:
	splitLine.erase(std::remove(splitLine.begin(), splitLine.end(), ""), splitLine.end());

        length = splitLine.size();

        elemType = stoi(splitLine[1]);
        // elemType == 1: 2-node   line element
        // elemType == 3: 4-node   quadrangle
        // elemType == 4: 4-node   tetrahedra

	int number_tags = stoi(splitLine[2]);

        if (elemType == 1)
        {
          int no1 = stoi(splitLine[3 + number_tags]) - 1;
          int no2 = stoi(splitLine[4 + number_tags]) - 1;

          if (method_type.compare("tlcpdi") == 0 ||
              method_type.compare("ulcpdi") == 0)
          {
	    xpc0[nc * l][0] = xpc[nc * l][0] = nodes[no1][0];
	    xpc0[nc * l][1] = xpc[nc * l][1] = nodes[no1][1];
	    xpc0[nc * l][2] = xpc[nc * l][2] = nodes[no1][2];

	    xpc0[nc * l + 1][0] = xpc[nc * l + 1][0] = nodes[no2][0];
	    xpc0[nc * l + 1][1] = xpc[nc * l + 1][1] = nodes[no2][1];
	    xpc0[nc * l + 1][2] = xpc[nc * l + 1][2] = nodes[no2][2];
	  }

	  x0[l][0] = x[l][0] = 0.5 * (nodes[no1][0] + nodes[no2][0]);
          x0[l][1] = x[l][1] = 0.5 * (nodes[no1][1] + nodes[no2][1]);
          x0[l][2] = x[l][2] = 0.5 * (nodes[no1][2] + nodes[no2][2]);
        }
        else if (elemType == 3)
        {

          // #ifdef DEBUG
          //           xcplot.clear();
          //           xcplot.resize(5, 0);
          //           ycplot.clear();
          //           ycplot.resize(5, 0);
          // #endif

          int no1 = stoi(splitLine[3 + number_tags]) - 1;
          int no2 = stoi(splitLine[4 + number_tags]) - 1;
          int no3 = stoi(splitLine[5 + number_tags]) - 1;
          int no4 = stoi(splitLine[6 + number_tags]) - 1;

          if (method_type.compare("tlcpdi") == 0 ||
              method_type.compare("ulcpdi") == 0)
          {
            xpc0[nc * l][0] = xpc[nc * l][0] = nodes[no1][0];
            xpc0[nc * l][1] = xpc[nc * l][1] = nodes[no1][1];
            xpc0[nc * l][2] = xpc[nc * l][2] = nodes[no1][2];

            xpc0[nc * l + 1][0] = xpc[nc * l + 1][0] = nodes[no2][0];
            xpc0[nc * l + 1][1] = xpc[nc * l + 1][1] = nodes[no2][1];
            xpc0[nc * l + 1][2] = xpc[nc * l + 1][2] = nodes[no2][2];

            xpc0[nc * l + 2][0] = xpc[nc * l + 2][0] = nodes[no3][0];
            xpc0[nc * l + 2][1] = xpc[nc * l + 2][1] = nodes[no3][1];
            xpc0[nc * l + 2][2] = xpc[nc * l + 2][2] = nodes[no3][2];

            xpc0[nc * l + 3][0] = xpc[nc * l + 3][0] = nodes[no4][0];
            xpc0[nc * l + 3][1] = xpc[nc * l + 3][1] = nodes[no4][1];
            xpc0[nc * l + 3][2] = xpc[nc * l + 3][2] = nodes[no4][2];
          }

          x0[l][0] = x[l][0] = 0.25 * (nodes[no1][0] + nodes[no2][0] +
                                         nodes[no3][0] + nodes[no4][0]);
          x0[l][1] = x[l][1] = 0.25 * (nodes[no1][1] + nodes[no2][1] +
                                         nodes[no3][1] + nodes[no4][1]);
          x0[l][2] = x[l][2] = 0.25 * (nodes[no1][2] + nodes[no2][2] +
                                         nodes[no3][2] + nodes[no4][2]);

          // vol0[l] = vol[l] = 0.5*(xpc[nc*l+0][0]*xpc[nc*l+1][1] -
          // xpc[nc*l+1][0]*xpc[nc*l+0][1]
          // 	    + xpc[nc*l+1][0]*xpc[nc*l+2][1] -
          // xpc[nc*l+2][0]*xpc[nc*l+1][1]
          // 	    + xpc[nc*l+2][0]*xpc[nc*l+3][1] -
          // xpc[nc*l+3][0]*xpc[nc*l+2][1]
          // 	    + xpc[nc*l+3][0]*xpc[nc*l+0][1] -
          // xpc[nc*l+0][0]*xpc[nc*l+3][1]);

          vol0[l] = vol[l] =
              0.5 *
              (nodes[no1][0] * nodes[no2][1] - nodes[no2][0] * nodes[no1][1] +
               nodes[no2][0] * nodes[no3][1] - nodes[no3][0] * nodes[no2][1] +
               nodes[no3][0] * nodes[no4][1] - nodes[no4][0] * nodes[no3][1] +
               nodes[no4][0] * nodes[no1][1] - nodes[no1][0] * nodes[no4][1]);
        }
        else if (elemType == 4)
        {

          int no1 = stoi(splitLine[3 + number_tags]) - 1;
          int no2 = stoi(splitLine[4 + number_tags]) - 1;
          int no3 = stoi(splitLine[5 + number_tags]) - 1;
          int no4 = stoi(splitLine[6 + number_tags]) - 1;

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

          x0[l][0] = x[l][0] = 0.25 * (x1 + x2 + x3 + x4);
          x0[l][1] = x[l][1] = 0.25 * (y1 + y2 + y3 + y4);
          x0[l][2] = x[l][2] = 0.25 * (z1 + z2 + z3 + z4);

          vol0[l] = vol[l] = (1 / 6) * (x21 * (y23 * z34 - y34 * z23) +
                                          x32 * (y34 * z12 - y12 * z34) +
                                          x43 * (y12 * z23 - y23 * z12));
        }
        else
        {
          cout << "Element type " << elemType << " not supported!!\n";
          error->one(FLERR, "");
        }

	// Check if the particle is inside the region:
	if (domain->inside_subdomain(x0[l][0], x0[l][1], x0[l][2])) {
	  l++;
	}
      }

      getline(file, line);
      if (line.compare("$EndElements") == 0)
      {
        if (universe->me == 0)
          cout << "Reading elements...done!\n";
        break;
      }
      else
	error->all(FLERR, "Unexpected line: " + line + ". $EndElements expected!!\n");
    }
  }

  if (np_local > l)
    grow(l);

  np_local = l; // Adjust np_local to account for the particles outside the domain

  if (universe->me == 0)
    cout << "np_local=" << np_local << endl;

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

  for (int i = 0; i < np_local; i++)
  {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    mbp[i].setZero();
    v_update[i].setZero();
    rho0[i] = rho[i] = mat->rho0;
    mass[i] = mat->rho0 * vol0[i];
    if (domain->axisymmetric == true) {
      mass[i] *= x0[i][0];
      vol0[i] *= x0[i][0];
    }
    vol[i] = vol0[i];

    eff_plastic_strain[i]      = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
    if (update->method->temp)
    {
      T[i] = T0;
      gamma[i] = 0;
      q[i].setZero();
    }
    ienergy[i] = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
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
  domain->np_local += np_local;
}

void Solid::compute_temperature_nodes(bool reset)
{
  double Ttemp;
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
    {
      grid->T[in] = 0;
    }

    if (grid->mass[in] > 0)
    {
      Ttemp = 0;

      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        Ttemp += wf_np[in][j] * mass[ip] * T[ip];
      }
      Ttemp /= grid->mass[in];
      grid->T[in] += Ttemp;
    }
  }
}

void Solid::compute_external_temperature_driving_forces_nodes(bool reset)
{
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    if (reset)
      grid->Qext[in] = 0;

    if (grid->mass[in] > 0)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        grid->Qext[in] += wf_np[in][j] * gamma[ip];
      }
    }
  }
}

void Solid::compute_internal_temperature_driving_forces_nodes()
{
  int ip, nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in = 0; in < nn; in++)
  {
    grid->Qint[in] = 0;
    for (int j = 0; j < numneigh_np[in]; j++)
    {
      ip = neigh_np[in][j];
      grid->Qint[in] += wfd_np[in][j].dot(q[ip]);

      if (domain->axisymmetric == true)
      {
        error->one(FLERR, "Temperature and axisymmetric not yet supported.\n");
        // ftemp[0] -= vol0PK1[ip](2, 2) * wf_np[in][j] / x0[ip][0];
      }
    }
  }
}

void Solid::update_particle_temperature()
{
  int in;
  for (int ip = 0; ip < np_local; ip++)
  {
    T[ip] = 0;
    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      // T[ip] += wf_pn[ip][j] * (grid->T_update[in] - grid->T[in]);
      T[ip] += wf_pn[ip][j] * grid->T_update[in];
    }
  }
}

void Solid::update_heat_flux(bool doublemapping)
{
  int in;

  vector<double> *Tn;

  if (doublemapping)
    Tn = &grid->T;
  else
    Tn = &grid->T_update;

  if (is_TL)
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      q[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        q[ip] -= wfd_pn[ip][j] * (*Tn)[in];
      }
      q[ip] *= vol0[ip] * mat->invcp * mat->kappa;
    }
  }
  else
  {
    for (int ip = 0; ip < np_local; ip++)
    {
      q[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        q[ip] -= wfd_pn[ip][j] * (*Tn)[in];
      }
      q[ip] *= vol[ip] * mat->invcp * mat->kappa;
    }
  }
}

void Solid::write_restart(ofstream *of)
{
  // Write solid bounds:
  of->write(reinterpret_cast<const char *>(&solidlo[0]), 3 * sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidhi[0]), 3 * sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidsublo[0]), 3 * sizeof(double));
  of->write(reinterpret_cast<const char *>(&solidsubhi[0]), 3 * sizeof(double));

  // Write number of particles:
  of->write(reinterpret_cast<const char *>(&np), sizeof(bigint));
  of->write(reinterpret_cast<const char *>(&np_local), sizeof(int));
  of->write(reinterpret_cast<const char *>(&nc), sizeof(int));

  // Write material's info:
  int iMat = material->find_material(mat->id);
  of->write(reinterpret_cast<const char *>(&iMat), sizeof(int));

  // Write cellsize:
  of->write(reinterpret_cast<const char *>(&grid->cellsize), sizeof(double));

  // Write particle's attributes:
  // cout << x[0](0) << ", " << x[0](1) << ", " << x[0](2) << endl;
  for (int ip = 0; ip < np_local; ip++)
  {
    of->write(reinterpret_cast<const char *>(&ptag[ip]), sizeof(tagint));
    of->write(reinterpret_cast<const char *>(&x0[ip]), sizeof(Eigen::Vector3d));
    of->write(reinterpret_cast<const char *>(&x[ip]), sizeof(Eigen::Vector3d));
    of->write(reinterpret_cast<const char *>(&v[ip]), sizeof(Eigen::Vector3d));
    of->write(reinterpret_cast<const char *>(&sigma[ip]), sizeof(Eigen::Matrix3d));
    of->write(reinterpret_cast<const char *>(&strain_el[ip]), sizeof(Eigen::Matrix3d));
    if (is_TL)
    {
      of->write(reinterpret_cast<const char *>(&vol0PK1[ip]), sizeof(Eigen::Matrix3d));
    }
    of->write(reinterpret_cast<const char *>(&F[ip]), sizeof(Eigen::Matrix3d));
    of->write(reinterpret_cast<const char *>(&J[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&vol0[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&rho0[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&eff_plastic_strain[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&eff_plastic_strain_rate[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&damage[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&damage_init[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&T[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&ienergy[ip]), sizeof(double));
    of->write(reinterpret_cast<const char *>(&mask[ip]), sizeof(int));
  }
}

void Solid::read_restart(ifstream *ifr)
{
  // Read solid bounds:
  ifr->read(reinterpret_cast<char *>(&solidlo[0]), 3 * sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidhi[0]), 3 * sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidsublo[0]), 3 * sizeof(double));
  ifr->read(reinterpret_cast<char *>(&solidsubhi[0]), 3 * sizeof(double));
  // cout << "solidlo=[" << solidlo[0] << "," << solidlo[1] << "," << solidlo[2] << endl;
  // cout << "solidhi=[" << solidhi[0] << "," << solidhi[1] << "," << solidhi[2] << endl;
  // cout << "solidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2] << endl;
  // cout << "solidsubhi=[" << solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2] << endl;

  // init();
  //  Read number of particles:
  ifr->read(reinterpret_cast<char *>(&np), sizeof(bigint));
  ifr->read(reinterpret_cast<char *>(&np_local), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&nc), sizeof(int));

  // Read material's info:
  int iMat = -1;
  ifr->read(reinterpret_cast<char *>(&iMat), sizeof(int));
  mat = &material->materials[iMat];

  // Read cellsize:
  ifr->read(reinterpret_cast<char *>(&grid->cellsize), sizeof(double));
  if (is_TL)
  {
    grid->init(solidlo, solidhi);
  }

  // Read particle's attributes:
  grow(np_local);

  for (int ip = 0; ip < np_local; ip++)
  {
    ifr->read(reinterpret_cast<char *>(&ptag[ip]), sizeof(tagint));
    ifr->read(reinterpret_cast<char *>(&x0[ip]), sizeof(Eigen::Vector3d));
    ifr->read(reinterpret_cast<char *>(&x[ip]), sizeof(Eigen::Vector3d));
    ifr->read(reinterpret_cast<char *>(&v[ip]), sizeof(Eigen::Vector3d));
    v_update[ip].setZero();
    a[ip].setZero();
    mbp[ip].setZero();
    f[ip].setZero();
    ifr->read(reinterpret_cast<char *>(&sigma[ip]), sizeof(Eigen::Matrix3d));
    ifr->read(reinterpret_cast<char *>(&strain_el[ip]), sizeof(Eigen::Matrix3d));
    if (is_TL)
    {
      ifr->read(reinterpret_cast<char *>(&vol0PK1[ip]), sizeof(Eigen::Matrix3d));
    }
    L[ip].setZero();
    ifr->read(reinterpret_cast<char *>(&F[ip]), sizeof(Eigen::Matrix3d));
    R[ip].setZero();
    D[ip].setZero();
    Finv[ip].setZero();
    Fdot[ip].setZero();
    ifr->read(reinterpret_cast<char *>(&J[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&vol0[ip]), sizeof(double));
    vol[ip] = J[ip] * vol0[ip];
    ifr->read(reinterpret_cast<char *>(&rho0[ip]), sizeof(double));
    rho[ip] = rho0[ip] / J[ip];
    mass[ip] = rho0[ip] * vol0[ip];
    ifr->read(reinterpret_cast<char *>(&eff_plastic_strain[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&eff_plastic_strain_rate[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&damage[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&damage_init[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&T[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&ienergy[ip]), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&mask[ip]), sizeof(int));
  }
  // cout << x[0](0) << ", " << x[0](1) << ", " << x[0](2) << endl;
}

void Solid::distribute_particles_by_domain()
{
  if (universe->nprocs == 1)
    return;

  int i, locator;
  // cout << "proc " << universe->me << " tagrange dst dom " << ptag[0] << " " << np_local << " " << ptag[np_local - 1] << endl;
  // evaluate the range of tags each proc uses to later determine the owner by particle tag
  tagrange.resize(universe->nprocs);
  tagint maxtag = ptag[np_local - 1] + 1;
  MPI_Allgather(&maxtag, 1, MPI_DOUBLE, tagrange.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

  // holds the data of particles to send for each process
  vector<vector<double>> send_prtcl(universe->nprocs);
  vector<int> sendcounts(universe->nprocs, 0);
  vector<int> displacements(universe->nprocs, 0);
  for (i = 0; i < np_local; ++i)
  {
    locator = domain->which_CPU_owns_me(x[i](0), x[i](1), x[i](2));
    if (locator != universe->me)
    {
      pack_particle(i, send_prtcl[locator]);
      // counts number of sent doubles
      sendcounts[locator] += comm_n;
    }
  }

  // calculate displacements
  for (i = 1; i < universe->nprocs; ++i)
    displacements[i] = displacements[i - 1] + sendcounts[i - 1];
  // stack the single vectors to one big sendbuf
  vector<double> sendbuf;
  for (auto v : send_prtcl)
    sendbuf.insert(sendbuf.begin(), v.begin(), v.end());

  // every process scatters its particles to every process in which domain the particle is located
  int recv_count, recv_p_count, iprt, index;
  vector<double> recvbuf;
  np_ghost = 0;
  for (i = 0; i < universe->nprocs; ++i)
  {
    if (i == universe->me)
    {
      MPI_Scatter(sendcounts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, i, MPI_COMM_WORLD);
      recvbuf.resize(recv_count);
      // cout << "proc " << universe->me << " iteration " << i << " recv count: "<<  recv_count << " buffsize: " << recvbuf.size() << endl;
      // for (int c = 0; c < sendcounts.size(); ++c)
      //   cout << "me: " << universe->me << " sendproc: " <<  i << " sendcounts: " << c << " my buff size: " << recvbuf.size() << "size of send_prtcl " << 0 << endl;
      MPI_Scatterv(sendbuf.data(), sendcounts.data(), displacements.data(), MPI_DOUBLE, recvbuf.data(), recv_count, MPI_DOUBLE, i, MPI_COMM_WORLD);
      // cout << "proc " << universe->me << " iteration " << i << " scatterv successfull!" << endl;
    }
    else
    {
      MPI_Scatter(NULL, 0, MPI_INT, &recv_count, 1, MPI_INT, i, MPI_COMM_WORLD);
      recvbuf.resize(recv_count);
      recv_p_count = recv_count / comm_n;
      np_ghost += recv_p_count;
      MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf.data(), recv_count, MPI_DOUBLE, i, MPI_COMM_WORLD);
      grow(np_local + np_ghost);
      // append new particles at the end of each vector
      for (iprt = 0; iprt < recv_p_count; ++iprt)
      {
        index = np_local + np_ghost - recv_p_count + iprt;
        // tags kommen erfolgreich an... werden aber falsch unpacked
        unpack_particle(index, iprt * comm_n, recvbuf);
      }
    }
    // a better implementation will replace every iterartion over np_local with np_local + np_ghost!
    np_local += np_ghost;
    // for (int j = np_local - np_ghost; j < np_local; ++j)
    //   cout << "proc " << universe->me << " ptag " << ptag[j] << endl;
    // cout << "proc " << universe->me << " distributed particles to domains! np_ghost: " << np_ghost << endl;
  }
}

void Solid::distribute_particles_by_process()
{
  if (universe->nprocs == 1)
    return;
  // vector<vector<double>> send_prtcl(universe->nprocs);
  vector<double> sendbuf;
  int i;
  tagint iproc;
  vector<int> send_counts(universe->nprocs, 0);
  vector<int> displacements(universe->nprocs, 0);
  // cout << "proc " << universe->me << " tagrange dst proc " << ptag[0] << " " << np_local << " " << ptag[np_local - 1] << endl;

  if (tagrange.size() != universe->nprocs)
    error->all(FLERR, "error: you need to call distribute by domain first.\n");

  // assuming that ghost particles keep the order of ptags
  for (i = np_local - np_ghost; i < np_local; ++i)
  {
    if (ptag[i] > tagrange[universe->nprocs - 1])
      error->all(FLERR, "Error: proc " + to_string(universe->me) + ": tag " + to_string(ptag[i]) + " out of tag range.");
    iproc = 0;
    // find process to which ptag belongs
    while (ptag[i] < tagrange[iproc] && iproc < universe->nprocs)
      ++iproc;
    send_counts[iproc] += comm_n;
    pack_particle(i, sendbuf);
  }

  // cout << "proc " << universe->me << " packed p by proc\n";
  for (i = 1; i < universe->nprocs; ++i)
    displacements[i] = displacements[i - 1] + send_counts[i - 1];
  // cout << "proc " << universe->me << " counts " << send_counts[i] << endl;

  // every particle is send to the process where it was created
  int recv_count, recv_p_count, start, ibuf;
  vector<double> recvbuf;
  for (i = 0; i < universe->nprocs; ++i)
  {
    if (i == universe->me)
    {
      MPI_Scatter(send_counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, i, MPI_COMM_WORLD);
      // cout << "proc " << universe->me << " scattered recv_count: " << recv_count << endl;
      recvbuf.resize(recv_count);
      MPI_Scatterv(sendbuf.data(), send_counts.data(), displacements.data(), MPI_DOUBLE, recvbuf.data(), recv_count, MPI_DOUBLE, i, MPI_COMM_WORLD);
      // cout << "proc " << universe->me << " scatterv successfull" << endl;
    }
    else
    {
      MPI_Scatter(NULL, 0, MPI_INT, &recv_count, 1, MPI_INT, i, MPI_COMM_WORLD);
      recvbuf.resize(recv_count);
      recv_p_count = recv_count / comm_n;
      MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuf.data(), recv_count, MPI_DOUBLE, i, MPI_COMM_WORLD);

      // now copy particles from buffer to matching ptag position
      vector<tagint>::iterator lower;
      tagint tag;
      // cout << "proc " << universe->me << " ordering tags." << endl;
      for (ibuf = 0; ibuf < recv_count; ibuf += comm_n)
      {
        tag = (tagint)recvbuf[ibuf];
        lower = lower_bound(ptag.begin(), ptag.end(), tag);
        // if (*lower != tag)
        // cout << "tag found: " << *lower << " tag searched: " << tag << endl;
        //   error->all(FLERR, "Error: could not find tag " + to_string(tag) + " in proc " + to_string(universe->me) + " \n");
        unpack_particle(lower - ptag.begin(), ibuf, recvbuf);
      }
      // cout << "proc " << universe->me << " copying particles in right order successfull. ibuff: " << ibuf << endl;
    }
  }
  // cout << "proc " << universe->me << " redistributed particles to procs\n";
  np_local -= np_ghost;
  np_ghost = 0;
  grow(np_local);
}

void Solid::get_particle_counts_and_displacements(int *counts, int *displacements, const int &root)
{
  // save counts and displacements of particles that are gathered to root
  int p, count_local;
  for (p = 0; p < universe->nprocs; ++p)
    counts[p] = displacements[p] = 0;
  count_local = np_local * comm_n;
  MPI_Gather(&count_local, 1, MPI_INT, counts, 1, MPI_INT, root, MPI_COMM_WORLD);
  // calc the particles of root by subtracting all other counts from np
  counts[root] = np * comm_n;
  for (p = 0; p < universe->nprocs; ++p)
  {
    if (p != root)
      counts[root] -= counts[p];
  }
  for (p = 1; p < universe->nprocs; ++p)
  {
    displacements[p] = displacements[p - 1] + counts[p - 1];
  }
}

void Solid::gather_particles(int &root)
{
  if (universe->nprocs == 1)
    return;
  auto counts_recv = new int[universe->nprocs];
  auto displacement = new int[universe->nprocs];

  get_particle_counts_and_displacements(counts_recv, displacement, root);

  vector<double> sendbuff, recvbuff;
  // store new particles in proc ROOT
  if (universe->me == root)
  {
    sendbuff.resize(np_local * comm_n);
    recvbuff.resize(np * comm_n);
    MPI_Gatherv(sendbuff.data(), np_local * comm_n, MPI_DOUBLE, recvbuff.data(),
                counts_recv, displacement, MPI_DOUBLE, root, MPI_COMM_WORLD);
    grow(np);

    // list from np_local to np to iterate over all particles except ROOT's
    vector<int> unpack_list(np - np_local);
    for (int i = 0; i < unpack_list.size(); ++i)
    {
      unpack_list[i] = (i + np_local) * comm_n;
    }
    unpack_particle(np_local, unpack_list, recvbuff);
  }
  else
  {
    for (int ip = 0; ip < np_local; ++ip)
    {
      pack_particle(ip, sendbuff);
    }
    MPI_Gatherv(sendbuff.data(), sendbuff.size(), MPI_DOUBLE, NULL,
                NULL, NULL, MPI_DOUBLE, root, MPI_COMM_WORLD);
  }

  delete counts_recv;
  delete displacement;
}

void Solid::scatter_particles(int &root)
{
  if (universe->nprocs == 1)
    return;
  auto counts_send = new int[universe->nprocs];
  auto displacement = new int[universe->nprocs];

  get_particle_counts_and_displacements(counts_send, displacement, root);
  vector<double> sendbuff, recvbuff;
  int buffsize;
  if (universe->me == root)
  {
    // pack all particles in buff
    for (int ip = 0; ip < np_local; ++ip)
    {
      pack_particle(ip, sendbuff);
    }
    np_local = counts_send[root] / comm_n;
    buffsize = np_local * comm_n;
    recvbuff.resize(buffsize);
    MPI_Scatterv(sendbuff.data(), counts_send, displacement, MPI_DOUBLE,
                 recvbuff.data(), buffsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
    grow(np_local);
  }
  else
  {
    buffsize = np_local * comm_n;
    recvbuff.resize(buffsize);
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, recvbuff.data(), buffsize, MPI_DOUBLE, root, MPI_COMM_WORLD);
    vector<int> unpack_list(np_local);
    int start = 0;
    for (int i = 0; i < np_local; ++i)
    {
      unpack_list[i] = i * comm_n;
    }
    // copy particles from buff to procs particles
    unpack_particle(start, unpack_list, recvbuff);
  }

  delete counts_send;
  delete displacement;
}

void Solid::surfmask_init(double cellsize)
{
  this->surf_cellsize = cellsize;
  bigint yinc, zinc, rsize;
  surfmask_dim[0] = (bigint)((domain->subhi[0] - domain->sublo[0]) / surf_cellsize);
  surfmask_dim[1] = (bigint)((domain->subhi[1] - domain->sublo[1]) / surf_cellsize);
  surfmask_dim[2] = (bigint)((domain->subhi[2] - domain->sublo[2]) / surf_cellsize);

  // if (domain->subhi[0] - domain->sublo[0] - surfmask_dim[0] * cellsize == 0)
  //   surfmask_dim[0]++;
  // if (domain->subhi[1] - domain->sublo[1] - surfmask_dim[1] * cellsize == 0)
  //   surfmask_dim[1]++;
  // if (domain->subhi[2] - domain->sublo[2] - surfmask_dim[2] * cellsize == 0)
  //   surfmask_dim[2]++;

  rsize = surfmask_dim[0] * surfmask_dim[1] * surfmask_dim[2];
  yinc = surfmask_dim[0];
  zinc = surfmask_dim[0] * surfmask_dim[1];

  surfmask = vector<int>(rsize, 0);
  p_in_cell = vector<vector<int>>(rsize);

  if (domain->dimension == 2) // 2D
  {
    // top, left, right, bottom cell
    surfmask_off[0] = +0 - yinc;
    surfmask_off[1] = -1;
    surfmask_off[2] = +1;
    surfmask_off[3] = +0 + yinc;
#if NEIGH2D == 6
    // diagonal cells
    surfmask_off[4] = -1 - yinc;
    surfmask_off[5] = +1 - yinc;
    surfmask_off[6] = -1 + yinc;
    surfmask_off[7] = +1 + yinc;
#endif
  }
  else if (domain->dimension == 3) // 3D
  {
    // face center cells
    surfmask_off[0] = +0 - zinc;
    surfmask_off[1] = +0 - yinc;
    surfmask_off[2] = -1;
    surfmask_off[3] = +1;
    surfmask_off[4] = +0 + yinc;
    surfmask_off[5] = +0 + zinc;

#if NEIGH3D >= 18
    // diagonal cells
    surfmask_off[6] = +0 - yinc - zinc;
    surfmask_off[7] = -1 - zinc;
    surfmask_off[8] = +1 - zinc;
    surfmask_off[9] = +0 + yinc - zinc;

    surfmask_off[10] = -1 - yinc;
    surfmask_off[11] = +1 - yinc;
    surfmask_off[12] = -1 + yinc;
    surfmask_off[13] = +1 + yinc;

    surfmask_off[14] = +0 - yinc + zinc;
    surfmask_off[15] = -1 + zinc;
    surfmask_off[16] = +1 + zinc;
    surfmask_off[17] = +0 + yinc + zinc;

#if NEIGH3D == 26
    // corner cells
    surfmask_off[18] = -1 - yinc - zinc;
    surfmask_off[19] = +1 - yinc - zinc;
    surfmask_off[20] = -1 + yinc - zinc;
    surfmask_off[21] = +1 + yinc - zinc;

    surfmask_off[22] = -1 - yinc + zinc;
    surfmask_off[23] = +1 - yinc + zinc;
    surfmask_off[24] = -1 + yinc + zinc;
    surfmask_off[25] = +1 + yinc + zinc;
#endif //! NEIGH3D >= 18
#endif //! NEIGH3D == 26
  }
}

void Solid::compute_surface_particles()
{
  if (np_local == 0)
    return;

  fill(is_surf.begin(), is_surf.end(), 0);
  fill(surfmask.begin(), surfmask.end(), 0);
  for (auto c : p_in_cell)
    c.clear();

  // cout << "proc " << universe->me << " cleared vectors" << endl;
  bigint i, rsize, coord, yinc, zinc;
  rsize = surfmask_dim[0] * surfmask_dim[1] * surfmask_dim[2];
  yinc = surfmask_dim[0];
  zinc = surfmask_dim[0] * surfmask_dim[1];
  // cout << "proc " << universe->me << " incs: 1 " << yinc << " "<<zinc << endl;
  // mark all positions where a particle is inside a cell on the grid
  bigint ix, iy, iz;
  for (i = 0; i < np_local; ++i)
  {
    if (domain->inside_subdomain(x[i](0), x[i](1), x[i](2)))
    {
      ix = (bigint)((x[i](0) - domain->sublo[0]) / surf_cellsize);
      if (ix == surfmask_dim[0])  // 
        --ix;
      iy = (bigint)((x[i](1) - domain->sublo[1]) / surf_cellsize);
      if (iy == surfmask_dim[1])
        --iy;
      iz = (bigint)((x[i](2) - domain->sublo[2]) / surf_cellsize);
      if (iz == surfmask_dim[2])
        --iz;
      coord = ix + iy * yinc + iz * zinc;
      surfmask[coord] = 1 << IS_SURF | 1 << HAS_PARTICLE;
      p_in_cell[coord].push_back(i);
    }
  }

  // mark grids as "no surface" if it has nneigh neighbours
  // neighbouring offsets and numbers
  int nneighmax, n_neigh, nreduced, flag;
  if (domain->dimension == 2)
    nneighmax = NEIGH2D;
  else if (domain->dimension == 3)
    nneighmax = NEIGH3D;

  do
  {
    nreduced = 0;
    for (i = 0; i < rsize; ++i)
    {
      if (surfmask[i] & 1 << IS_SURF)
      {
        n_neigh = 0;
        iz = i / zinc;
        iy = (i % zinc) / yinc;
        ix = i % yinc;
        // is cell not on border?
        flag = (ix < 1) << XLO_BOUND;
        flag |= (iy < 1) << YLO_BOUND;
        flag |= (iz < 1 && domain->dimension == 3) << ZLO_BOUND;
        flag |= (ix >= surfmask_dim[0] - 1) << XHI_BOUND;
        flag |= (iy >= surfmask_dim[1] - 1) << YHI_BOUND;
        flag |= ((iz >= surfmask_dim[2] - 1) && (domain->dimension == 3)) << ZHI_BOUND;

        // switch (flag & IS_BOUND)
        // {
        // case C_LLL:
        //   surfmask[i] = 1;
        //   break;
        // case C_HLL:
        //   surfmask[i] = 2;
        //   break;
        // case C_LHL:
        //   surfmask[i] = 3;
        //   break;
        // case C_HHL:
        //   surfmask[i] = 4;
        //   break;
        // case C_LLH:
        //   surfmask[i] = 5;
        //   break;
        // case C_HLH:
        //   surfmask[i] = 6;
        //   break;
        // case C_LHH:
        //   surfmask[i] = 7;
        //   break;
        // case C_HHH:
        //   surfmask[i] = 8;
        //   break;

        // case E_XLO_YLO:
        //   surfmask[i] = 9;
        //   break;
        // case E_XHI_YLO:
        //   surfmask[i] = 10;
        //   break;
        // case E_XLO_YHI:
        //   surfmask[i] = 11;
        //   break;
        // case E_XHI_YHI:
        //   surfmask[i] = 12;
        //   break;

        // case E_XLO_ZLO:
        //   surfmask[i] = 13;
        //   break;
        // case E_XHI_ZLO:
        //   surfmask[i] = 14;
        //   break;
        // case E_XLO_ZHI:
        //   surfmask[i] = 15;
        //   break;
        // case E_XHI_ZHI:
        //   surfmask[i] = 16;
        //   break;

        // case E_YLO_ZLO:
        //   surfmask[i] = 17;
        //   break;
        // case E_YHI_ZLO:
        //   surfmask[i] = 18;
        //   break;
        // case E_YLO_ZHI:
        //   surfmask[i] = 19;
        //   break;
        // case E_YHI_ZHI:
        //   surfmask[i] = 20;
        //   break;
        // case 1 << XLO_BOUND:
        //   surfmask[i] = 21;
        //   break;
        // case 1 << XHI_BOUND:
        //   surfmask[i] = 22;
        //   break;
        // case 1 << YLO_BOUND:
        //   surfmask[i] = 23;
        //   break;
        // case 1 << YHI_BOUND:
        //   surfmask[i] = 24;
        //   break;
        // case 1 << ZLO_BOUND:
        //   surfmask[i] = 25;
        //   break;
        // case 1 << ZHI_BOUND:
        //   surfmask[i] = 26;
        //   break;
        // default:
        //   surfmask[i] = 0;
        //   break;
        // }

        if (flag & IS_BOUND) // particle is on bound
        {
          surfmask[i] |= flag; // surfmask now signals if this is a bound particle
        }

        if (surfmask[i] & IS_SURF)
        {
          for (int neigh = 0; neigh < nneighmax; ++neigh)
          {
            if (p_in_cell[i + surfmask_off[neigh]].empty() == false)
            {
              ++n_neigh;
            }
          }
          if (n_neigh >= nneighmax)
          {
            // is surf bit is set to zero
            surfmask[i] &= ~(1 << IS_SURF);
            ++nreduced;
          }
        }
      }
    }
  } while (nreduced > 0);

  for (i = 0; i < rsize; ++i)
  {
    if (surfmask[i] & 1 << HAS_PARTICLE)
    {
      for (auto pid : p_in_cell[i])
      {
        is_surf[pid] = surfmask[i];
      }
    }
  }
}

// switch (flag)
// {
// case C_LLL:
//   break;
// case C_HLL:
//   break;
// case C_LHL:
//   break;
// case C_HHL:
//   break;
// case C_LLH:
//   break;
// case C_HLH:
//   break;
// case C_LHH:
//   break;
// case C_HHH:
//   break;

// case E_XLO_YLO:
//   break;
// case E_XHI_YLO:
//   break;
// case E_XLO_YHI:
//   break;
// case E_XHI_YHI:
//   break;

// case E_XLO_ZLO:
//   break;
// case E_XHI_ZLO:
//   break;
// case E_XLO_ZHI:
//   break;
// case E_XHI_ZHI:
//   break;

// case E_YLO_ZLO:
//   break;
// case E_YHI_ZLO:
//   break;
// case E_YLO_ZHI:
//   break;
// case E_YHI_ZHI:
//   break;

// default:
//   for (int neigh = 0; neigh < nneighmax; ++neigh)
//   {
//     if (p_in_cell[i + surfmask_off[neigh]].empty() == false)
//     {
//       ++n_neigh;
//     }
//   }
//   if (n_neigh == nneighmax)
//   {
//     // cout << "found surrounded particle: " << i << endl;
//     surfmask[i] = 0;
//     ++nreduced;
//   }
//   break;
// }