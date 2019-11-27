#include "solid.h"
#include "domain.h"
#include "input.h"
#include "material.h"
#include "memory.h"
#include "method.h"
#include "mpm.h"
#include "mpm_math.h"
#include "update.h"
#include "var.h"
#include <Eigen/Eigen>
#include <math.h>
#include <omp.h>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

#define SQRT_3_OVER_2 1.224744871 // sqrt(3.0/2.0)

Solid::Solid(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == NULL)
  {
    cout << "Error: a method should be defined before creating a solid!"
         << endl;
    exit(1);
  }

  if (args.size() < 2)
  {
    cout << "Error: solid command not enough arguments. " << endl;
    for (auto &x : usage)
      cout << x.second;
    exit(1);
  }

  if (usage.find(args[1]) == usage.end())
  {
    cout << "Error, keyword \033[1;31m" << args[1] << "\033[0m unknown!\n";
    for (auto &x : usage)
      cout << x.second;
    exit(1);
  }

  if (args.size() < Nargs.find(args[1])->second)
  {
    cout << "Error: not enough arguments.\n";
    cout << usage.find(args[1])->second;
    exit(1);
  }

  cout << "Creating new solid with ID: " << args[0] << endl;

  method_type = update->method_type;
  id          = args[0];

  np = 0;

  x = x0 = NULL;
  rp = rp0 = NULL;
  xpc = xpc0 = NULL;

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {
    nc = pow(2, domain->dimension);
  }
  else
    nc = 0;

  v = v_update = NULL;

  a = NULL;

  sigma = strain_el = vol0PK1 = L = F = R = U = D = Finv = Fdot = Di = NULL;

  mb = f = NULL;

  J   = NULL;
  vol = vol0 = NULL;
  rho = rho0              = NULL;
  mass                    = NULL;
  eff_plastic_strain      = NULL;
  eff_plastic_strain_rate = NULL;
  damage                  = NULL;
  damage_init             = NULL;
  T                       = NULL;
  ienergy                 = NULL;
  mask                    = NULL;

  mat = NULL;

  if (method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    grid = new Grid(mpm);

  else
    grid = domain->grid;

  numneigh_pn = numneigh_np = NULL;

  neigh_pn = neigh_np = NULL;

  wf_pn = wf_np = wf_pn_corners = NULL;
  wfd_pn = wfd_np = NULL;

  dtCFL = 1.0e22;
  vtot  = 0;

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

#ifdef DEBUG
  plt::axis("equal");
  plt::save("debug.png");
  plt::close();
  // exit(1);
#endif
}

Solid::~Solid()
{
  if (x0 != NULL)
    delete x0;
  if (x != NULL)
    delete x;
  if (rp0 != NULL)
    delete rp0;
  if (rp != NULL)
    delete rp;
  if (xpc0 != NULL)
    delete xpc0;
  if (xpc != NULL)
    delete xpc;
  if (v != NULL)
    delete v;
  if (v_update != NULL)
    delete v_update;
  if (a != NULL)
    delete a;
  if (mb != NULL)
    delete mb;
  if (f != NULL)
    delete f;
  if (sigma != NULL)
    delete sigma;
  if (strain_el != NULL)
    delete strain_el;
  if (vol0PK1 != NULL)
    delete vol0PK1;
  if (L != NULL)
    delete L;
  if (F != NULL)
    delete F;
  if (R != NULL)
    delete R;
  if (U != NULL)
    delete U;
  if (D != NULL)
    delete D;
  if (Finv != NULL)
    delete Finv;
  if (Fdot != NULL)
    delete Fdot;
  if (Di != NULL)
    delete Di;

  memory->destroy(J);
  memory->destroy(vol);
  memory->destroy(vol0);
  memory->destroy(rho);
  memory->destroy(rho0);
  memory->destroy(mass);
  memory->destroy(eff_plastic_strain);
  memory->destroy(eff_plastic_strain_rate);
  memory->destroy(damage);
  memory->destroy(damage_init);
  memory->destroy(T);
  memory->destroy(ienergy);
  memory->destroy(mask);

  if (method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    delete grid;

  delete[] numneigh_pn;
  delete[] numneigh_np;

  delete[] neigh_pn;
  delete[] neigh_np;

  delete[] wf_pn;
  delete[] wf_np;
  if (wf_pn_corners != NULL)
    delete[] wf_pn_corners;

  delete[] wfd_pn;
  delete[] wfd_np;
}

void Solid::init()
{
  cout << "Bounds for " << id << ":\n";
  cout << "xlo xhi: " << solidlo[0] << " " << solidhi[0] << endl;
  cout << "ylo yhi: " << solidlo[1] << " " << solidhi[1] << endl;
  cout << "zlo zhi: " << solidlo[2] << " " << solidhi[2] << endl;

  // Calculate total volume and mass
  vtot = 0.;

  double mtot(0.);

  for (int ip = 0; ip < np; ip++)
  {
    vtot += vol[ip];
    mtot += mass[ip];
  }

  cout << "Solid " << id << " total volume = " << vtot << endl;
  cout << "Solid " << id << " total mass = " << mtot << endl;

  if (grid->nnodes == 0)
    grid->init(solidlo, solidhi);

  if (np == 0)
  {
    cout << "Error: solid does not have any particles" << endl;
    exit(1);
  }
  else
  {
    bigint nnodes = grid->nnodes;

    numneigh_pn = new int[np]();
    neigh_pn    = new vector<int>[np];
    wf_pn       = new vector<double>[np];
    if (nc != 0)
      wf_pn_corners = new vector<double>[nc * np];
    wfd_pn = new vector<Vector3d>[np];

    if (nnodes)
    {
      numneigh_np = new int[nnodes]();
      neigh_np    = new vector<int>[nnodes];
      wf_np       = new vector<double>[nnodes];
      wfd_np      = new vector<Vector3d>[nnodes];
    }
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
  if (args->end() < it + 3)
  {
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it)
  {
    int iMat = material->find_material(*it);

    if (iMat == -1)
    {
      cout << "Error: could not find material named " << *it << endl;
      exit(1);
    }

    mat = &material->materials[iMat]; // point mat to the right material

    it++;

    if (grid->cellsize == 0)
      grid->setup(*it); // set the grid cellsize

    it++;
    T0 = input->parsev(*it); // set initial temperature

    it++;

    if (it != args->end())
    {
      cout << "Error: too many arguments" << endl;
      for (auto &x : usage)
        cout << x.second;
      exit(1);
    }
  }
}

void Solid::grow(int nparticles)
{
  np = nparticles;

  string str;

  str = "solid-" + id + ":x0";
  cout << "Growing " << str << endl;
  if (x0 == NULL)
    x0 = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: x0 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":x";
  cout << "Growing " << str << endl;
  if (x == NULL)
    x = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: x already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {

    if (update->method->style == 0)
    { // CPDI-R4
      str = "solid-" + id + ":rp0";
      cout << "Growing " << str << endl;
      if (rp0 == NULL)
        rp0 = new Eigen::Vector3d[domain->dimension * np];
      else
      {
        cout << "Error: rp0 already exists, I don't know how to grow it!\n";
        exit(1);
      }

      str = "solid-" + id + ":rp";
      cout << "Growing " << str << endl;
      if (rp == NULL)
        rp = new Eigen::Vector3d[domain->dimension * np];
      else
      {
        cout << "Error: rp already exists, I don't know how to grow it!\n";
        exit(1);
      }
    }
    if (update->method->style == 1)
    { // CPDI-Q4
      str = "solid-" + id + ":xpc0";
      cout << "Growing " << str << endl;
      if (xpc0 == NULL)
        xpc0 = new Eigen::Vector3d[nc * np];
      else
      {
        cout << "Error: xpc0 already exists, I don't know how to grow it!\n";
        exit(1);
      }
      for (int i = 0; i < nc * np; i++)
        xpc0[i].setZero();

      str = "solid-" + id + ":xpc";
      cout << "Growing " << str << endl;
      if (xpc == NULL)
        xpc = new Eigen::Vector3d[nc * np];
      else
      {
        cout << "Error: xpc already exists, I don't know how to grow it!\n";
        exit(1);
      }
      for (int i = 0; i < nc * np; i++)
        xpc[i].setZero();
    }
  }

  str = "solid-" + id + ":v";
  cout << "Growing " << str << endl;
  if (v == NULL)
    v = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: v already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":v_update";
  cout << "Growing " << str << endl;
  if (v_update == NULL)
    v_update = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: v_update already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":a";
  cout << "Growing " << str << endl;
  if (a == NULL)
    a = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: a already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":mb";
  cout << "Growing " << str << endl;
  if (mb == NULL)
    mb = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: mb already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":f";
  cout << "Growing " << str << endl;
  if (f == NULL)
    f = new Eigen::Vector3d[np];
  else
  {
    cout << "Error: f already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (sigma == NULL)
    sigma = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: sigma already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (strain_el == NULL)
    strain_el = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: strain_el already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (vol0PK1 == NULL)
    vol0PK1 = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: vol0PK1 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (L == NULL)
    L = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: L already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (F == NULL)
    F = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: F already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (R == NULL)
    R = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: R already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (U == NULL)
    U = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: U already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (D == NULL)
    D = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: D already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (Finv == NULL)
    Finv = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: Finv already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (Fdot == NULL)
    Fdot = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: Fdot already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (Di == NULL)
    Di = new Eigen::Matrix3d[np];
  else
  {
    cout << "Error: Di already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":vol0";
  cout << "Growing " << str << endl;
  vol0 = memory->grow(vol0, np, str);

  str = "solid-" + id + ":vol";
  cout << "Growing " << str << endl;
  vol = memory->grow(vol, np, str);

  str = "solid-" + id + ":rho0";
  cout << "Growing " << str << endl;
  rho0 = memory->grow(rho0, np, str);

  str = "solid-" + id + ":rho";
  cout << "Growing " << str << endl;
  rho = memory->grow(rho, np, str);

  str = "solid-" + id + ":mass";
  cout << "Growing " << str << endl;
  mass = memory->grow(mass, np, str);

  str = "solid-" + id + ":eff_plastic_strain";
  cout << "Growing " << str << endl;
  eff_plastic_strain = memory->grow(eff_plastic_strain, np, str);

  str = "solid-" + id + ":eff_plastic_strain_rate";
  cout << "Growing " << str << endl;
  eff_plastic_strain_rate = memory->grow(eff_plastic_strain_rate, np, str);

  str = "solid-" + id + ":damage";
  cout << "Growing " << str << endl;
  damage = memory->grow(damage, np, str);

  str = "solid-" + id + ":damage_init";
  cout << "Growing " << str << endl;
  damage_init = memory->grow(damage_init, np, str);

  str = "solid-" + id + ":T";
  cout << "Growing " << str << endl;
  T = memory->grow(T, np, str);

  str = "solid-" + id + ":ienergy";
  cout << "Growing " << str << endl;
  ienergy = memory->grow(ienergy, np, str);

  str = "solid-" + id + ":mask";
  cout << "Growing " << str << endl;
  mask = memory->grow(mask, np, str);

  for (int i = 0; i < np; i++)
    mask[i] = 1;

  str = "solid-" + id + ":J";
  cout << "Growing " << str << endl;
  J = memory->grow(J, np, str);
}

void Solid::compute_mass_nodes(bool reset)
{
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (reset)
      grid->mass[in] = 0;

    // if (grid->rigid[in] && !mat->rigid) continue;

    for (int j = 0; j < numneigh_np[in]; j++)
    {
      ip = neigh_np[in][j];
      grid->mass[in] += wf_np[in][j] * mass[ip];
      // if (in==0) {
      // cout << "compute_mass_nodes:\ttag=" << in << "\tptag = " << ip <<
      // "\tmass[ip]=" << mass[ip] << "\tphi=" << wf_np[in][j] << "\tmassn=" <<
      // grid->mass[in] << endl;
      // }
    }
  }
  return;
}

void Solid::compute_velocity_nodes(bool reset)
{
  Eigen::Vector3d *vn        = grid->v;
  Eigen::Vector3d *vn_update = grid->v_update;
  Eigen::Vector3d vtemp, vtemp_rigid;
  double *massn = grid->mass;
  double mass_rigid;
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (reset)
    {
      vn[in].setZero();
      vn_update[in].setZero();
    }

    if (mat->rigid)
      mass_rigid = 0;

    if (massn[in] > 0)
    {
      vtemp.setZero();
      if (mat->rigid)
        vtemp_rigid.setZero();

      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        vtemp += (wf_np[in][j] * mass[ip]) * v[ip];

        if (grid->rigid[in] && mat->rigid)
        {
          vtemp_rigid += wf_np[in][j] * v[ip];
          mass_rigid += wf_np[in][j];
        }
        // vn[in] += (wf_np[in][j] * mass[ip]) * v[ip]/ massn[in];
      }
      vtemp /= massn[in];
      vn[in] += vtemp;

      if (mat->rigid && mass_rigid > 1.0e-12)
        vn_update[in] += vtemp_rigid / mass_rigid;
      if (isnan(vn_update[in](0)))
        cout << "in=" << in << "\tvn=[" << vn[in][0] << ", " << vn[in][1]
             << ", " << vn[in][2] << "]\tvp=[" << v[ip][0] << ", " << v[ip][1]
             << ", " << v[ip][2] << "],\tvn_update=[" << vn_update[in][0]
             << ", " << vn_update[in][1] << ", " << vn_update[in][2] << "]\n";
    }
  }
}

void Solid::compute_velocity_nodes_APIC(bool reset)
{
  Eigen::Vector3d *x0n = grid->x0;
  Eigen::Vector3d *vn  = grid->v;
  double *massn        = grid->mass;
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (reset)
      vn[in].setZero();

    if (grid->rigid[in] && !mat->rigid)
      continue;

    if (massn[in] > 0)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        vn[in] += (wf_np[in][j] * mass[ip]) *
                  (v[ip] + Fdot[ip] * (x0n[in] - x0[ip])) / massn[in];
      }
    }
  }
}

void Solid::compute_external_forces_nodes(bool reset)
{
  Eigen::Vector3d *mbn = grid->mb;
  double *massn        = grid->mass;
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (reset)
      mbn[in].setZero();

    if (grid->rigid[in])
      continue;

    if (massn[in] > 0)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        mbn[in] += wf_np[in][j] * mb[ip];
      }
    }
  }
}

void Solid::compute_internal_forces_nodes_TL()
{
  Eigen::Vector3d *fn = grid->f;
  Eigen::Vector3d ftemp;
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (grid->rigid[in])
    {
      fn[in].setZero();
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

    fn[in] = ftemp;
  }
}

void Solid::compute_internal_forces_nodes_UL(bool reset)
{
  Eigen::Vector3d *fn = grid->f;
  double *massn       = grid->mass;
  int ip;

  for (int in = 0; in < grid->nnodes; in++)
  {
    if (reset)
      fn[in].setZero();

    if (grid->rigid[in])
      continue;

    for (int j = 0; j < numneigh_np[in]; j++)
    {
      ip = neigh_np[in][j];
      fn[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
    }

    if (domain->axisymmetric == true)
    {
      for (int j = 0; j < numneigh_np[in]; j++)
      {
        ip = neigh_np[in][j];
        fn[in][0] -= vol[ip] * (sigma[ip](2, 2) * wf_np[in][j] / x[ip][0]);
      }
    }
  }
}

void Solid::compute_particle_velocities_and_positions()
{
  Eigen::Vector3d *vn_update = grid->v_update;

  vector<Eigen::Vector3d> vc_update;
  vc_update.resize(nc);

  int in;

  bool update_corners, ul;

  if (update->method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    ul = false;
  else
    ul = true;

  if ((method_type.compare("tlcpdi") == 0 ||
       method_type.compare("ulcpdi") == 0) &&
      (update->method->style == 1))
  {
    update_corners = true;
  }
  else
    update_corners = false;

  for (int ip = 0; ip < np; ip++)
  {

    v_update[ip].setZero();
    if (update_corners)
      for (int i = 0; i < nc; i++)
        vc_update[i].setZero();

    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * vn_update[in];
      x[ip] += update->dt * wf_pn[ip][j] * vn_update[in];
      if (isnan(x[ip](0)))
        cout << "ip=" << ip << "\tx=[" << x[ip](0) << "," << x[ip](1) << ","
             << x[ip](2) << "]\tin=" << in << "\tvn_update=["
             << vn_update[in](0) << "," << vn_update[in](1) << ","
             << vn_update[in](2) << "]\twf_pn=" << wf_pn[ip][j] << "\n";

      if (update_corners)
      {
        for (int ic = 0; ic < nc; ic++)
        {
          vc_update[ic] += wf_pn_corners[nc * ip + ic][j] * vn_update[in];
        }
      }
    }

    if (ul)
    {
      // Check if the particle is within the box's domain:
      if (domain->inside(x[ip]) == 0)
      {
        cout << "Error: Particle " << ip << " left the domain ("
             << domain->boxlo[0] << "," << domain->boxhi[0] << ","
             << domain->boxlo[1] << "," << domain->boxhi[1] << ","
             << domain->boxlo[2] << "," << domain->boxhi[2] << ",):\n"
             << x[ip] << endl;
        exit(1);
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

  Eigen::Vector3d *vn_update = grid->v_update;
  Eigen::Vector3d *vn        = grid->v;

  int in;

  for (int ip = 0; ip < np; ip++)
  {
    a[ip].setZero();
    if (mat->rigid)
      continue;
    for (int j = 0; j < numneigh_pn[ip]; j++)
    {
      in = neigh_pn[ip][j];
      a[ip] += wf_pn[ip][j] * (vn_update[in] - vn[in]);
    }
    a[ip] *= inv_dt;
    // if (ip==234)
    //   cout << "ip=" << ip << "\ta=[" << a[ip](0) << "," << a[ip](1) << "," <<
    //   a[ip](2) << "]\n";
    f[ip] = a[ip] / mass[ip];
  }
}

void Solid::update_particle_velocities(double FLIP)
{

  if (mat->rigid)
    return;

  for (int ip = 0; ip < np; ip++)
  {
    v[ip] = (1 - FLIP) * v_update[ip] + FLIP * (v[ip] + update->dt * a[ip]);
  }
}

void Solid::compute_rate_deformation_gradient_TL()
{
  if (mat->rigid)
    return;

  int in;
  Eigen::Vector3d *vn = grid->v;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
      }
    }
  }
  else if (domain->dimension == 2)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];

        if (domain->axisymmetric == true)
          Fdot[ip](2, 2) += vn[in][0] * wf_pn[ip][j] / x0[ip][0];
      }
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        Fdot[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        Fdot[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        Fdot[ip](0, 2) += vn[in][0] * wfd_pn[ip][j][2];
        Fdot[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        Fdot[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
        Fdot[ip](1, 2) += vn[in][1] * wfd_pn[ip][j][2];
        Fdot[ip](2, 0) += vn[in][2] * wfd_pn[ip][j][0];
        Fdot[ip](2, 1) += vn[in][2] * wfd_pn[ip][j][1];
        Fdot[ip](2, 2) += vn[in][2] * wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_MUSL()
{
  if (mat->rigid)
    return;

  int in;
  Eigen::Vector3d *vn = grid->v;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
      }
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == false))
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        L[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
      }
    }
  }
  else if ((domain->dimension == 2) && (domain->axisymmetric == true))
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        L[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
        L[ip](2, 2) += vn[in][0] * wf_pn[ip][j] / x[ip][0];
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
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        L[ip](0, 2) += vn[in][0] * wfd_pn[ip][j][2];
        L[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
        L[ip](1, 2) += vn[in][1] * wfd_pn[ip][j][2];
        L[ip](2, 0) += vn[in][2] * wfd_pn[ip][j][0];
        L[ip](2, 1) += vn[in][2] * wfd_pn[ip][j][1];
        L[ip](2, 2) += vn[in][2] * wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_USL()
{
  if (mat->rigid)
    return;

  int in;
  Eigen::Vector3d *vn = grid->v_update;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
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
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        L[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
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
        L[ip](0, 0) += vn[in][0] * wfd_pn[ip][j][0];
        L[ip](0, 1) += vn[in][0] * wfd_pn[ip][j][1];
        L[ip](0, 2) += vn[in][0] * wfd_pn[ip][j][2];
        L[ip](1, 0) += vn[in][1] * wfd_pn[ip][j][0];
        L[ip](1, 1) += vn[in][1] * wfd_pn[ip][j][1];
        L[ip](1, 2) += vn[in][1] * wfd_pn[ip][j][2];
        L[ip](2, 0) += vn[in][2] * wfd_pn[ip][j][0];
        L[ip](2, 1) += vn[in][2] * wfd_pn[ip][j][1];
        L[ip](2, 2) += vn[in][2] * wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_deformation_gradient()
{
  if (mat->rigid)
    return;

  int in;
  Eigen::Vector3d *xn  = grid->x;
  Eigen::Vector3d *x0n = grid->x0;
  Eigen::Vector3d dx;
  Eigen::Matrix3d Ftemp, eye;
  eye.setIdentity();

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = xn[in] - x0n[in];
        Ftemp(0, 0) += dx[0] * wfd_pn[ip][j][0];
      }
      F[ip](0, 0) = Ftemp(0, 0) + 1;
    }
  }
  else if (domain->dimension == 2)
  {
    for (int ip = 0; ip < np; ip++)
    {
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = xn[in] - x0n[in];
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
    for (int ip = 0; ip < np; ip++)
    {
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = xn[in] - x0n[in];
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
  Eigen::Vector3d *x0n = grid->x0;
  // Eigen::Vector3d *vn = grid->v;
  Eigen::Vector3d *vn = grid->v_update;
  Eigen::Vector3d dx;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = x0n[in] - x0[ip];
        Fdot[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  }
  else if (domain->dimension == 2)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = x0n[in] - x0[ip];
        Fdot[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
        Fdot[ip](0, 1) += vn[in][0] * dx[1] * wf_pn[ip][j];
        Fdot[ip](1, 0) += vn[in][1] * dx[0] * wf_pn[ip][j];
        Fdot[ip](1, 1) += vn[in][1] * dx[1] * wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip = 0; ip < np; ip++)
    {
      Fdot[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = x0n[in] - x0[ip];
        Fdot[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
        Fdot[ip](0, 1) += vn[in][0] * dx[1] * wf_pn[ip][j];
        Fdot[ip](0, 2) += vn[in][0] * dx[2] * wf_pn[ip][j];
        Fdot[ip](1, 0) += vn[in][1] * dx[0] * wf_pn[ip][j];
        Fdot[ip](1, 1) += vn[in][1] * dx[1] * wf_pn[ip][j];
        Fdot[ip](1, 2) += vn[in][1] * dx[2] * wf_pn[ip][j];
        Fdot[ip](2, 0) += vn[in][2] * dx[0] * wf_pn[ip][j];
        Fdot[ip](2, 1) += vn[in][2] * dx[1] * wf_pn[ip][j];
        Fdot[ip](2, 2) += vn[in][2] * dx[2] * wf_pn[ip][j];
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
  Eigen::Vector3d *x0n = grid->x0;
  // Eigen::Vector3d *vn = grid->v;
  Eigen::Vector3d *vn = grid->v_update;
  Eigen::Vector3d dx;

  if (domain->dimension == 1)
  {
    for (int ip = 0; ip < np; ip++)
    {
      L[ip].setZero();
      for (int j = 0; j < numneigh_pn[ip]; j++)
      {
        in = neigh_pn[ip][j];
        dx = x0n[in] - x0[ip];
        L[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
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
        dx = x0n[in] - x0[ip];
        L[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
        L[ip](0, 1) += vn[in][0] * dx[1] * wf_pn[ip][j];
        L[ip](1, 0) += vn[in][1] * dx[0] * wf_pn[ip][j];
        L[ip](1, 1) += vn[in][1] * dx[1] * wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
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
        dx = x0n[in] - x0[ip];
        L[ip](0, 0) += vn[in][0] * dx[0] * wf_pn[ip][j];
        L[ip](0, 1) += vn[in][0] * dx[1] * wf_pn[ip][j];
        L[ip](0, 2) += vn[in][0] * dx[2] * wf_pn[ip][j];
        L[ip](1, 0) += vn[in][1] * dx[0] * wf_pn[ip][j];
        L[ip](1, 1) += vn[in][1] * dx[1] * wf_pn[ip][j];
        L[ip](1, 2) += vn[in][1] * dx[2] * wf_pn[ip][j];
        L[ip](2, 0) += vn[in][2] * dx[0] * wf_pn[ip][j];
        L[ip](2, 1) += vn[in][2] * dx[1] * wf_pn[ip][j];
        L[ip](2, 2) += vn[in][2] * dx[2] * wf_pn[ip][j];
      }
      L[ip] *= Di[ip];
    }
  }
}

void Solid::update_deformation_gradient()
{
  if (mat->rigid)
    return;

  bool status, tl, lin, nh, vol_cpdi;
  Eigen::Matrix3d U;
  Eigen::Matrix3d eye;
  eye.setIdentity();

  if (update->method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    tl = true;
  else
    tl = false;

  // if (mat->type == material->constitutive_model::LINEAR) lin = true;
  // else lin = false;

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

  for (int ip = 0; ip < np; ip++)
  {

    if (tl)
      F[ip] += update->dt * Fdot[ip];
    else
      F[ip] = (eye + update->dt * L[ip]) * F[ip];

    Finv[ip] = F[ip].inverse();

    if (J[ip] < 0.0)
    {
      cout << "Error: J[" << ip << "]<0.0 == " << J[ip] << endl;
      cout << "F[" << ip << "]:" << endl << F[ip] << endl;
      exit(1);
    }

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
    rho[ip] = rho0[ip] / J[ip];

    if (!nh)
    {
      // Only done if not Neo-Hookean:

      // TLMPM. L is computed from Fdot
      if (tl)
      {
        // cout << Finv[ip] << endl;
        // pocout << Fdot[ip] << endl;
        L[ip] = Fdot[ip] * Finv[ip];

        // if  ( domain->axisymmetric == true ) L[ip](2, 2) += vn[in][0] *
        // wf_pn[ip][j] / x[ip][0];
      }
      // else
      //   Fdot[ip] = L[ip]*F[ip];

      status = PolDec(
          F[ip], R[ip], U,
          false); // polar decomposition of the deformation gradient, F = R * U

      if (tl)
        D[ip] = 0.5 * (R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip]);
      else
        D[ip] = 0.5 * (L[ip] + L[ip].transpose());

      if (!status)
      {
        cout << "Polar decomposition of deformation gradient failed for "
                "particle "
             << ip << ".\n";
        cout << "F:" << endl << F[ip] << endl;
        cout << "timestep" << endl << update->ntimestep << endl;
        exit(1);
      }
    }

    // strain_increment[ip] = update->dt * D[ip];
  }
}

void Solid::update_stress()
{
  if (mat->rigid)
    return;

  max_p_wave_speed = 0;
  double pH, plastic_strain_increment, flow_stress;
  Matrix3d eye, sigma_dev, FinvT, PK1, strain_increment;
  bool lin, tl, nh, fluid, temp;

  if (mat->type == material->constitutive_model::LINEAR)
    lin = true;
  else
    lin = false;

  if (mat->type == material->constitutive_model::NEO_HOOKEAN)
    nh = true;
  else
    nh = false;

  if (update->method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    tl = true;
  else
    tl = false;

  if (mat->temp.size())
    temp = true;
  else
    temp = false;

  eye.setIdentity();
  plastic_strain_increment = 0;
  sigma_dev.setZero();

  //# pragma omp parallel for
  for (int ip = 0; ip < np; ip++)
  {
    if (lin)
    {
      strain_increment = update->dt * D[ip];
      strain_el[ip] += strain_increment;
      sigma[ip] += 2 * mat->G * strain_increment +
                   mat->lambda * strain_increment.trace() * eye;
      if (temp)
        sigma_dev = Deviator(sigma[ip]);
      if (tl)
        vol0PK1[ip] = vol0[ip] * J[ip] *
                      (R[ip] * sigma[ip] * R[ip].transpose()) *
                      Finv[ip].transpose();
    }
    else if (nh)
    {

      // Neo-Hookean material:
      FinvT       = Finv[ip].transpose();
      PK1         = mat->G * (F[ip] - FinvT) + mat->lambda * log(J[ip]) * FinvT;
      vol0PK1[ip] = vol0[ip] * PK1;
      sigma[ip]   = 1.0 / J[ip] * (F[ip] * PK1.transpose());
      if (temp)
        sigma_dev = Deviator(sigma[ip]);
      strain_el[ip] =
          0.5 * (F[ip].transpose() * F[ip] - eye); // update->dt * D[ip];
    }
    else
    {

      mat->eos->compute_pressure(pH, ienergy[ip], J[ip], rho[ip], T[ip],
                                 damage[ip]);
      sigma_dev = mat->strength->update_deviatoric_stress(
          sigma[ip], D[ip], plastic_strain_increment, eff_plastic_strain[ip],
          eff_plastic_strain_rate[ip], damage[ip], T[ip]);

      eff_plastic_strain[ip] += plastic_strain_increment;

      // // compute a characteristic time over which to average the plastic
      // strain
      double tav = 1000 * grid->cellsize / mat->signal_velocity;
      eff_plastic_strain_rate[ip] -=
          eff_plastic_strain_rate[ip] * update->dt / tav;
      eff_plastic_strain_rate[ip] += plastic_strain_increment / tav;
      eff_plastic_strain_rate[ip] = MAX(0.0, eff_plastic_strain_rate[ip]);

      if (mat->damage != NULL)
        mat->damage->compute_damage(damage_init[ip], damage[ip], pH, sigma_dev,
                                    eff_plastic_strain_rate[ip],
                                    plastic_strain_increment, T[ip]);
      sigma[ip] = -pH * eye + sigma_dev;

      if (damage[ip] > 1e-10)
      {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye +
            sigma_dev / (mat->G * (1 - damage[ip]));
      }
      else
      {
        strain_el[ip] =
            (update->dt * D[ip].trace() + strain_el[ip].trace()) / 3.0 * eye;
      }

      if (tl)
      {
        vol0PK1[ip] = vol0[ip] * J[ip] *
                      (R[ip] * sigma[ip] * R[ip].transpose()) *
                      Finv[ip].transpose();
      }
    }

    if (temp)
    {
      flow_stress = SQRT_3_OVER_2 * sigma_dev.norm();
      for (int itemp = 0; itemp < mat->temp.size(); itemp++)
        mat->temp[itemp]->compute_temperature(T[ip], flow_stress,
                                              plastic_strain_increment);
    }
  }

  double min_h_ratio = 1.0e22;
  double four_third  = 1.333333333333333333333333333333333333333;
  for (int ip = 0; ip < np; ip++)
  {
    max_p_wave_speed =
        MAX(max_p_wave_speed,
            sqrt((mat->K + four_third * mat->G) / rho[ip]) +
                MAX(MAX(abs(v[ip](0)), abs(v[ip](1))), abs(v[ip](2))));

    min_h_ratio =
        MIN(min_h_ratio, F[ip](0, 0) * F[ip](0, 0) + F[ip](0, 1) * F[ip](0, 1) +
                             F[ip](0, 2) * F[ip](0, 2));
    min_h_ratio =
        MIN(min_h_ratio, F[ip](1, 0) * F[ip](1, 0) + F[ip](1, 1) * F[ip](1, 1) +
                             F[ip](1, 2) * F[ip](1, 2));
    min_h_ratio =
        MIN(min_h_ratio, F[ip](2, 0) * F[ip](2, 0) + F[ip](2, 1) * F[ip](2, 1) +
                             F[ip](2, 2) * F[ip](2, 2));

    if (std::isnan(max_p_wave_speed))
    {
      cout << "Error: max_p_wave_speed is nan with ip=" << ip
           << ", rho[ip]=" << rho[ip] << ", K=" << mat->K << ", G=" << mat->G
           << endl;
      exit(1);
    }
    else if (max_p_wave_speed < 0.0)
    {
      cout << "Error: max_p_wave_speed= " << max_p_wave_speed
           << " with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K
           << ", G=" << mat->G << endl;
      exit(1);
    }
  }

  dtCFL = MIN(dtCFL, grid->cellsize * sqrt(min_h_ratio) / max_p_wave_speed);
  if (std::isnan(dtCFL))
  {
    cout << "Error: dtCFL = " << dtCFL << "\n";
    cout << "max_p_wave_speed = " << max_p_wave_speed
         << ", grid->cellsize=" << grid->cellsize << endl;
    exit(1);
  }
}

void Solid::compute_inertia_tensor(string form_function)
{

  // int in;
  // Eigen::Vector3d *x0n = grid->x0;
  // Eigen::Vector3d *vn_update = grid->v_update;
  // Eigen::Vector3d dx;

  // if (domain->dimension == 2) {
  //   for (int ip=0; ip<np; ip++){
  //     Di[ip].setZero();
  //     for (int j=0; j<numneigh_pn[ip]; j++){
  // 	in = neigh_pn[ip][j];
  // 	dx = x0n[in] - x0[ip];
  // 	Di[ip](0,0) += wf_pn[ip][j]*(dx[0]*dx[0]);
  // 	Di[ip](0,1) += wf_pn[ip][j]*(dx[0]*dx[1]);
  // 	Di[ip](1,1) += wf_pn[ip][j]*(dx[1]*dx[1]);
  //     }
  //     Di[ip](1,0) = Di[ip](0,1);
  //     cout << "Di[" << ip << "]=\n" << Di[ip] << endl;
  //   }
  // } else if (domain->dimension == 3) {
  //   for (int ip=0; ip<np; ip++){
  //     Di[ip].setZero();
  //     for (int j=0; j<numneigh_pn[ip]; j++){
  // 	in = neigh_pn[ip][j];
  // 	dx = x0n[in] - x0[ip];
  // 	Di[ip](0,0) += wf_pn[ip][j]*(dx[0]*dx[0]);
  // 	Di[ip](0,1) += wf_pn[ip][j]*(dx[0]*dx[1]);
  // 	Di[ip](0,2) += wf_pn[ip][j]*(dx[0]*dx[2]);
  // 	Di[ip](1,1) += wf_pn[ip][j]*(dx[1]*dx[1]);
  // 	Di[ip](1,2) += wf_pn[ip][j]*(dx[1]*dx[2]);
  // 	Di[ip](2,2) += wf_pn[ip][j]*(dx[2]*dx[2]);
  //     }
  //     Di[ip](1,0) = Di[ip](0,1);
  //     Di[ip](2,1) = Di[ip](1,2);
  //     Di[ip](2,0) = Di[ip](0,2);
  //     cout << "Di[" << ip << "]=\n" << Di[ip] << endl;
  //   }
  // }

  // Di[ip] = Di[ip].inverse();

  Eigen::Matrix3d eye;
  eye.setIdentity();

  double cellsizeSqInv = 1.0 / (grid->cellsize * grid->cellsize);

  for (int ip = 0; ip < np; ip++)
  {
    if (form_function.compare("linear") == 0)
    {
      // If the form function is linear:
      Di[ip] = 16.0 / 3.0 * cellsizeSqInv * eye;
    }
    else if (form_function.compare("quadratic-spline") == 0)
    {
      // If the form function is a quadratic spline:
      Di[ip] = 4.0 * cellsizeSqInv * eye;
    }
    else if (form_function.compare("cubic-spline") == 0)
    {
      // If the form function is a cubic spline:
      Di[ip] = 3.0 * cellsizeSqInv * eye;
    }
    else if (form_function.compare("Bernstein-quadratic") == 0)
      Di[ip] = 12.0 * cellsizeSqInv * eye;
    // cout << "Di[" << ip << "]=\n" << Di[ip] << endl;
  }
}

void Solid::copy_particle(int i, int j)
{
  x0[j]                      = x0[i];
  x[j]                       = x[i];
  v[j]                       = v[i];
  v_update[j]                = v[i];
  a[j]                       = a[i];
  mb[j]                      = mb[i];
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
  T[j]                       = T[i];
  ienergy[j]                 = ienergy[i];
  sigma[j]                   = sigma[i];
  vol0PK1[j]                 = vol0PK1[i];
  L[j]                       = L[i];
  F[j]                       = F[i];
  R[j]                       = R[i];
  U[j]                       = U[i];
  D[j]                       = D[i];
  Finv[j]                    = Finv[i];
  Fdot[j]                    = Fdot[i];
  J[j]                       = J[i];

  if (method_type.compare("tlcpdi") == 0 || method_type.compare("ulcpdi") == 0)
  {
    if (update->method->style == 0)
    { // CPDI-R4
      for (int id = 0; id < domain->dimension; id++)
      {
        rp0[domain->dimension * j + id] = rp0[domain->dimension * i + id];
        rp[domain->dimension * j + id]  = rp[domain->dimension * i + id];
      }
    }

    if (update->method->style == 1)
    { // CPDI-Q4
      for (int ic = 0; ic < nc; ic++)
      {
        xpc0[nc * j + ic] = xpc0[nc * i + ic];
        xpc[nc * j + ic]  = xpc[nc * i + ic];
      }
    }
  }
}

void Solid::populate(vector<string> args)
{

  cout << "Solid delimitated by region ID: " << args[2] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1)
  {
    cout << "Error: region ID " << args[2] << " does not exist" << endl;
    exit(1);
  }

  vector<double> limits = domain->regions[iregion]->limits();

  solidlo[0] = limits[0];
  solidhi[0] = limits[1];
  solidlo[1] = limits[2];
  solidhi[1] = limits[3];
  solidlo[2] = limits[4];
  solidhi[2] = limits[5];

  // Calculate total number of particles np:
  int nx, ny, nz;
  double delta;
  double hdelta;
  double Lx, Ly, Lz;

  double *boundlo;

  delta = grid->cellsize;

  if (grid->nnodes == 0)
  {
    // The grid will be ajusted to the solid's domain (good for TLMPM),
    // so all particles created will lie in the region:

    boundlo = solidlo;

    // and we need to create the corresponding grid:
    grid->init(solidlo, solidhi);

    Lx = solidhi[0] - solidlo[0];
    if (domain->dimension >= 2)
      Ly = solidhi[1] - solidlo[1];
    if (domain->dimension == 3)
      Lz = solidhi[2] - solidlo[2];
  }
  else
  {
    // The grid is most likely bigger than the solid's domain (good for ULMPM),
    // so all particles created won't lie in the region, they will need to be
    // checked:

    boundlo = domain->boxlo;

    Lx = domain->boxhi[0] - domain->boxlo[0];
    if (domain->dimension >= 2)
      Ly = domain->boxhi[1] - domain->boxlo[1];
    if (domain->dimension == 3)
      Lz = domain->boxhi[2] - domain->boxlo[2];
  }

  nx = (int)(Lx / delta);
  while (nx * delta <= Lx - 0.5 * delta)
    nx++;

  if (domain->dimension >= 2)
  {
    ny = (int)(Ly / delta);
    while (ny * delta <= Ly - 0.5 * delta)
      ny++;
  }
  else
  {
    ny = 1;
  }

  if (domain->dimension == 3)
  {
    nz = (int)(Lz / delta);
    while (nz * delta <= Lz - 0.5 * delta)
      nz++;
  }
  else
  {
    nz = 1;
  }

  np = nx * ny * nz;

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

  int np_per_cell = (int)input->parsev(args[3]);

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

    if (domain->dimension == 1)
      nip = 2;
    else if (domain->dimension == 2)
      nip = 4;
    else
      nip = 8;

    if (nc == 0)
      xi = 0.5 / sqrt(3.0);
    else
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
    if (domain->dimension == 1)
      nip = 3;
    else if (domain->dimension == 2)
      nip = 9;
    else
      nip = 27;

    intpoints = {-xi, -xi, -xi, -xi, 0,   -xi, -xi, xi,  -xi, 0,   -xi, -xi,
                 0,   0,   -xi, 0,   xi,  -xi, xi,  -xi, -xi, xi,  0,   -xi,
                 xi,  xi,  -xi, -xi, -xi, 0,   -xi, 0,   0,   -xi, xi,  0,
                 0,   -xi, 0,   0,   0,   0,   0,   xi,  0,   xi,  -xi, 0,
                 xi,  0,   0,   xi,  xi,  0,   -xi, -xi, xi,  -xi, 0,   xi,
                 -xi, xi,  xi,  0,   -xi, xi,  0,   0,   xi,  0,   xi,  xi,
                 xi,  -xi, xi,  xi,  0,   xi,  xi,  xi,  xi};
  }
  else
  {
    cout << "Error: solid command 4th argument should be 1,  2 or 3, but "
         << (int)input->parsev(args[3]) << "received.\n";
    exit(1);
  }

  np *= nip;
  mass_ /= (double)nip;
  vol_ /= (double)nip;

  // Allocate the space in the vectors for np particles:
  grow(np);

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

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        for (int ip = 0; ip < nip; ip++)
        {

          if (l >= np)
          {
            cout << "Error in Solid::populate(), exceeding the allocated "
                    "number of particles.\n";
            cout << "l = " << l << endl;
            cout << "np = " << np << endl;
            exit(1);
          }

          x0[l][0] = x[l][0] =
              boundlo[0] + delta * (i + 0.5 + intpoints[3 * ip + 0]);
          x0[l][1] = x[l][1] =
              boundlo[1] + delta * (j + 0.5 + intpoints[3 * ip + 1]);
          if (dim == 3)
            x0[l][2] = x[l][2] =
                boundlo[2] + delta * (k + 0.5 + intpoints[3 * ip + 2]);
          else
            x0[l][2] = x[l][2] = 0;

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
          // Check if the particle is inside the region:
          if (domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2]) ==
              1)
            l++;
        }
      }
    }
  }

  np = l; // Adjust np to account for the particles outside the domain
  cout << "np=" << np << endl;

  for (int i = 0; i < np; i++)
  {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    mb[i].setZero();
    v_update[i].setZero();
    vol0[i] = vol[i] = vol_;
    rho0[i] = rho[i] = mat->rho0;

    if (domain->axisymmetric == true)
    {
      mass[i] = mass_ * x0[i][0];
      vol0[i] = vol[i] = mass[i] / rho0[i];
    }
    else
      mass[i] = mass_;

    eff_plastic_strain[i]      = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i]                  = 0;
    damage_init[i]             = 0;
    T[i]                       = T0;
    ienergy[i]                 = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    U[i].setZero();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
    Di[i].setZero();

    J[i] = 1;
  }

  if (l != np)
  {
    cout << "Error l=" << l << " != np=" << np << endl;
    exit(1);
  }
}

void Solid::update_particle_domain()
{
  int dim = domain->dimension;

  if (update->method->style == 0)
  { // CPDI-R4
    for (int ip = 0; ip < np; ip++)
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

#ifdef DEBUG
  std::vector<double> x2plot, y2plot;
  std::vector<double> xcplot, ycplot;
#endif

  ifstream file(fileName, std::ios::in);

  if (!file)
  {
    cout << "Error: unable to open mesh file " << fileName << endl;
    exit(1);
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
        exit(1);
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
          exit(1);
        }

        if (domain->dimension == 2 && xn[2] != 0.)
        {
          cout << "Error: node " << id << " has non zero z component.\n";
          exit(1);
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

#ifdef DEBUG
          xcplot.clear();
          xcplot.resize(5, 0);
          ycplot.clear();
          ycplot.resize(5, 0);
#endif

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

#ifdef DEBUG
          xcplot[0] = xpc0[nc * ie][0];
          ycplot[0] = xpc0[nc * ie][1];
          xcplot[1] = xpc0[nc * ie + 1][0];
          ycplot[1] = xpc0[nc * ie + 1][1];
          xcplot[2] = xpc0[nc * ie + 2][0];
          ycplot[2] = xpc0[nc * ie + 2][1];
          xcplot[3] = xpc0[nc * ie + 3][0];
          ycplot[3] = xpc0[nc * ie + 3][1];
          xcplot[4] = xpc0[nc * ie][0];
          ycplot[4] = xpc0[nc * ie][1];
          plt::plot(xcplot, ycplot, "r-");
#endif

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

#ifdef DEBUG
          x2plot.push_back(x0[ie][0]);
          y2plot.push_back(x0[ie][1]);
#endif
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
          exit(1);
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
    mb[i].setZero();
    v_update[i].setZero();
    rho0[i] = rho[i]           = mat->rho0;
    mass[i]                    = mat->rho0 * vol0[i];
    eff_plastic_strain[i]      = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i]                  = 0;
    damage_init[i]             = 0;
    T[i]                       = T0;
    ienergy[i]                 = 0;
    strain_el[i].setZero();
    sigma[i].setZero();
    vol0PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    U[i].setZero();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
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

#ifdef DEBUG
  plt::plot(x2plot, y2plot, ".");
#endif
}
