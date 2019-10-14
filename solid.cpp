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

#include <mpi.h>
#include "mpm.h"
#include "solid.h"
#include "memory.h"
#include "update.h"
#include "method.h"
#include "domain.h"
#include "universe.h"
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include "mpm_math.h"
#include <math.h>
#include "input.h"
#include "var.h"
#include <omp.h>
#include "error.h"

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == NULL) {
    error->all(FLERR, "Error: a method should be defined before creating a solid!\n");
  }
  if (args.size() < 3) {
    error->all(FLERR, "Error: solid command not enough arguments.\n");
  }

  cout << "Creating new solid with ID: " << args[0] << endl;

  method_style = update->method_style;
  id = args[0];

  np = 0;

  if (update->method->is_CPDI) {
    nc = pow(2, domain->dimension);
  } else nc = 0;

  mat = NULL;

  if (update->method->is_TL) grid = new Grid(mpm);
  else grid = domain->grid;

  numneigh_pn = numneigh_np = NULL;

  neigh_pn = neigh_np = NULL;

  wf_pn = wf_np = NULL;
  wfd_pn = wfd_np = NULL;

  dtCFL = 1.0e22;
  vtot = 0;

  // Set material and cellsize:
  options(&args, args.begin()+3);

  // Create particles:
  populate(args);
}

Solid::~Solid()
{
  if (method_style.compare("tlmpm") == 0) delete grid;

  delete [] numneigh_pn;
  delete [] numneigh_np;

  delete [] neigh_pn;
  delete [] neigh_np;

  delete [] wf_pn;
  delete [] wf_np;

  delete [] wfd_pn;
  delete [] wfd_np;
}


void Solid::init()
{
  cout << "Bounds for " << id << ":\n";
  cout << "xlo xhi: " << solidlo[0] << " " << solidhi[0] << endl;
  cout << "ylo yhi: " << solidlo[1] << " " << solidhi[1] << endl;
  cout << "zlo zhi: " << solidlo[2] << " " << solidhi[2] << endl;

  // Calculate total volume:
  double vtot_local = 0;
  for (int ip=0; ip<np_local; ip++){
    vtot_local += vol[ip];
  }

  MPI_Allreduce(&vtot_local, &vtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

  cout << "Solid " << id << " total volume = " << vtot << endl;

  if (grid->nnodes_local == 0) grid->init(solidlo, solidhi);

  if (np_local == 0) {
    error->one(FLERR,"Error: solid does not have any particles.\n");
  } else {
      bigint nnodes = grid->nnodes_local + grid->nnodes_ghost;

      numneigh_pn = new int[np_local]();
      neigh_pn = new vector<int>[np_local];
      wf_pn = new vector<double>[np_local];
      wfd_pn = new vector< Vector3d >[np_local];

      if (nnodes) {
	numneigh_np = new int[nnodes]();
	neigh_np = new vector<int>[nnodes];
	wf_np = new vector<double>[nnodes];
	wfd_np = new vector< Vector3d >[nnodes];
      }
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
  if (args->end() < it+2) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it) {
    int iMat = material->find_material(*it);

    if (iMat == -1){
      cout << "Error: could not find material named " << *it << endl;
      error->all(FLERR,"\n");
    }

    mat = &material->materials[iMat]; // point mat to the right material

    it++;

    if (grid->cellsize == 0) grid->setup(*it); // set the grid cellsize

    it++;

    if (it != args->end()) {
      error->all(FLERR,"Error: too many arguments.\n");
    }
  }
}


void Solid::grow(int nparticles){
  //np_local = nparticles;

  string str;

  str = "solid-" + id + ":ptag";
  cout << "Growing " << str << endl;
  // if (ptag == NULL) ptag = new tagint[nparticles];
  ptag.resize(nparticles);

  str = "solid-" + id + ":x0";
  cout << "Growing " << str << endl;
  x0.resize(nparticles);

  str = "solid-" + id + ":x";
  cout << "Growing " << str << endl;
  x.resize(nparticles);

  if (method_style.compare("tlcpdi") == 0
      || method_style.compare("ulcpdi") == 0) {

    str = "solid-" + id + ":rp0";
    cout << "Growing " << str << endl;
    // if (rp0 == NULL) rp0 = new Eigen::Vector3d[domain->dimension*nparticles];
    rp0.resize(nparticles);

    str = "solid-" + id + ":rp";
    cout << "Growing " << str << endl;
    // if (rp == NULL) rp = new Eigen::Vector3d[domain->dimension*nparticles];
    rp.resize(nparticles);
  }

  if (method_style.compare("tlcpdi2") == 0
      || method_style.compare("ulcpdi2") == 0) {

    str = "solid-" + id + ":xpc0";
    cout << "Growing " << str << endl;
    // if (xpc0 == NULL) xpc0 = new Eigen::Vector3d[nc*nparticles];
    xpc0.resize(nparticles);

    str = "solid-" + id + ":xpc";
    cout << "Growing " << str << endl;
    // if (xpc == NULL) xpc = new Eigen::Vector3d[nc*nparticles];
    xpc.resize(nparticles);
  }

  str = "solid-" + id + ":v";
  cout << "Growing " << str << endl;
  // if (v == NULL) v = new Eigen::Vector3d[nparticles];
  v.resize(nparticles);

  str = "solid-" + id + ":v_update";
  cout << "Growing " << str << endl;
  // if (v_update == NULL) v_update = new Eigen::Vector3d[nparticles];
  v_update.resize(nparticles);

  str = "solid-" + id + ":a";
  cout << "Growing " << str << endl;
  // if (a == NULL) a = new Eigen::Vector3d[nparticles];
  a.resize(nparticles);

  str = "solid-" + id + ":mbp";
  cout << "Growing " << str << endl;
  mbp.resize(nparticles);

  str = "solid-" + id + ":f";
  cout << "Growing " << str << endl;
  // if (f == NULL) f = new Eigen::Vector3d[nparticles];
  f.resize(nparticles);

  // if (sigma == NULL) sigma = new Eigen::Matrix3d[nparticles];
  sigma.resize(nparticles);

  // if (strain_el == NULL) strain_el = new Eigen::Matrix3d[nparticles];
  strain_el.resize(nparticles);

  // if (vol0PK1 == NULL) vol0PK1 = new Eigen::Matrix3d[nparticles];
  vol0PK1.resize(nparticles);

  // if (L == NULL) L = new Eigen::Matrix3d[nparticles];
  L.resize(nparticles);

  // if (F == NULL) F = new Eigen::Matrix3d[nparticles];
  F.resize(nparticles);

  // if (R == NULL) R = new Eigen::Matrix3d[nparticles];
  R.resize(nparticles);

  // if (U == NULL) U = new Eigen::Matrix3d[nparticles];
  U.resize(nparticles);

  // if (D == NULL) D = new Eigen::Matrix3d[nparticles];
  D.resize(nparticles);

  // if (Finv == NULL) Finv = new Eigen::Matrix3d[nparticles];
  Finv.resize(nparticles);

  // if (Fdot == NULL) Fdot = new Eigen::Matrix3d[nparticles];
  Fdot.resize(nparticles);

  // if (Di == NULL) Di = new Eigen::Matrix3d[nparticles];
  Di.resize(nparticles);


  str = "solid-" + id + ":vol0";
  cout << "Growing " << str << endl;
  vol0.resize(nparticles);

  str = "solid-" + id + ":vol";
  cout << "Growing " << str << endl;
  vol.resize(nparticles);

  str = "solid-" + id + ":rho0";
  cout << "Growing " << str << endl;
  rho0.resize(nparticles);

  str = "solid-" + id + ":rho";
  cout << "Growing " << str << endl;
  rho.resize(nparticles);

  str = "solid-" + id + ":mass";
  cout << "Growing " << str << endl;
  mass.resize(nparticles);

  str = "solid-" + id + ":eff_plastic_strain";
  cout << "Growing " << str << endl;
  eff_plastic_strain.resize(nparticles);

  str = "solid-" + id + ":eff_plastic_strain_rate";
  cout << "Growing " << str << endl;
  eff_plastic_strain_rate.resize(nparticles);

  str = "solid-" + id + ":damage";
  cout << "Growing " << str << endl;
  damage.resize(nparticles);

  str = "solid-" + id + ":damage_init";
  cout << "Growing " << str << endl;
  damage_init.resize(nparticles);

  str = "solid-" + id + ":mask";
  cout << "Growing " << str << endl;
  mask.resize(nparticles);

  str = "solid-" + id + ":J";
  cout << "Growing " << str << endl;
  J.resize(nparticles);
}

void Solid::compute_mass_nodes(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++){
    if (reset) grid->mass[in] = 0;

    for (int j=0; j<numneigh_np[in];j++){
      ip = neigh_np[in][j];
      grid->mass[in] += wf_np[in][j] * mass[ip];
      // if (grid->ntag[in]==5) {
      // 	cout << "compute_mass_nodes:\ttag=" << grid->ntag[in] << "\tptag = " << ptag[ip] << "\tmass[ip]=" << mass[ip] << "\tphi=" << wf_np[in][j] << "\tmassn=" << grid->mass[in] << endl;
      // }
    }
  }
  return;
}

void Solid::compute_velocity_nodes(bool reset)
{
  Eigen::Vector3d vtemp;
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++) {
    if (reset) grid->v[in].setZero();
    if (grid->mass[in] > 0) {
      vtemp.setZero();
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	vtemp += (wf_np[in][j] * mass[ip]) * v[ip];
	// if (grid->ntag[in]==31 || grid->ntag[in]==64)
	//   printf("proc %d - ntag=%d, ptag=%d, ip=%d, wf=%.10e, mp=%.10e, vyp=%.10e, vtemp=%.10e\n", universe->me, grid->ntag[in], ptag[ip], ip, wf_np[in][j], mass[ip], v[ip][1], vtemp[1]);
	//(*vn)[in] += (wf_np[in][j] * mass[ip]) * v[ip]/ massn[in];
      }
      vtemp /= grid->mass[in];
      grid->v[in] += vtemp;
    }
  }
}

void Solid::compute_velocity_nodes_APIC(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++) {
    if (reset) grid->v[in].setZero();
    if (grid->mass[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	grid->v[in] += (wf_np[in][j] * mass[ip]) * (v[ip] + Fdot[ip]*(grid->x0[in] - x0[ip]))/ grid->mass[in];
      }
    }
  }
}

void Solid::compute_external_forces_nodes(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++) {
    if (reset)  grid->mb[in].setZero();
    if ( grid->mass[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	 grid->mb[in] += wf_np[in][j] *  mbp[ip];
	 // if (grid->ntag[in]==5706)
	 //   printf("proc %d - ntag=%d, ptag=%d, ip=%d, wf=%.10e,\tmbp=[%.10e, %.10e, %.10e],\tmb=[%.10e, %.10e, %.10e]\n", universe->me, grid->ntag[in], ptag[ip], ip, wf_np[in][j], mbp[ip][0], mbp[ip][1], mbp[ip][2], grid->mb[in][0], grid->mb[in][1], grid->mb[in][2]);
      }
    }
  }
}

void Solid::compute_internal_forces_nodes_TL()
{
  Eigen::Vector3d ftemp;
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++) {
    ftemp.setZero();
    for (int j=0; j<numneigh_np[in];j++){
      ip = neigh_np[in][j];
      ftemp -= vol0PK1[ip] * wfd_np[in][j];
    }
    grid->f[in] = ftemp;
  }
}

void Solid::compute_internal_forces_nodes_UL(bool reset)
{
  int ip;
  int nn = grid->nnodes_local + grid->nnodes_ghost;

  for (int in=0; in<nn; in++) {
    if (reset) grid->f[in].setZero();
    for (int j=0; j<numneigh_np[in];j++){
      ip = neigh_np[in][j];
      grid->f[in] -= vol[ip] * (sigma[ip] * wfd_np[in][j]);
    }
  }
}


void Solid::compute_particle_velocities()
{
  int in;

  for (int ip=0; ip<np_local; ip++){
    v_update[ip].setZero();
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * grid->v_update[in];
      // if (ptag[ip]==80||ptag[ip]==85)
      // 	printf("proc %d - ntag=%d, ptag=%d, ip=%d, wf=%.10e, vyp=%.10e, vyn=%.10e\n", universe->me, grid->ntag[in], ptag[ip], ip, wf_pn[ip][j], v_update[ip][1], grid->v_update[in][1]);
    }
  }
}

void Solid::compute_particle_acceleration()
{
  double inv_dt = 1.0/update->dt;

  int in;

  for (int ip=0; ip<np_local; ip++){
    a[ip].setZero();
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      a[ip] += wf_pn[ip][j] * (grid->v_update[in] - grid->v[in]);
    }
    a[ip] *= inv_dt;
    f[ip] = a[ip] / mass[ip];
  }
}

void Solid::update_particle_position()
{
  bool ul;

  if (update->method_style.compare("tlmpm") != 0) ul = true;
  else ul = false;

  for (int ip=0; ip<np_local; ip++) {
    x[ip] += update->dt*v_update[ip];
    if (ul) {
      // Check if the particle is within the box's domain:
      if (domain->inside(x[ip]) == 0) {
	error->all(FLERR,"Error: Particle " + to_string(ip) + " left the domain (" +
		   to_string(domain->boxlo[0]) + ","+ to_string(domain->boxhi[0]) + "," +
		   to_string(domain->boxlo[1]) + ","+ to_string(domain->boxhi[1]) + "," +
		   to_string(domain->boxlo[2]) + ","+ to_string(domain->boxhi[2]) + ",):\n");
      }
    }
  }
}

void Solid::update_particle_velocities(double FLIP)
{
  for (int ip=0; ip<np_local; ip++) {
    v[ip] = (1 - FLIP) * v_update[ip] + FLIP*(v[ip] + update->dt*a[ip]);
  }
}

void Solid::compute_rate_deformation_gradient_TL()
{
  int in;
  vector<Eigen::Vector3d> *vn = &grid->v;

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	Fdot[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
      }
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	Fdot[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	Fdot[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	Fdot[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	Fdot[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	Fdot[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	Fdot[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	Fdot[ip](0,2) += (*vn)[in][0]*wfd_pn[ip][j][2];
	Fdot[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	Fdot[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
	Fdot[ip](1,2) += (*vn)[in][1]*wfd_pn[ip][j][2];
	Fdot[ip](2,0) += (*vn)[in][2]*wfd_pn[ip][j][0];
	Fdot[ip](2,1) += (*vn)[in][2]*wfd_pn[ip][j][1];
	Fdot[ip](2,2) += (*vn)[in][2]*wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::compute_rate_deformation_gradient_UL_MUSL()
{
  int in;
  vector<Eigen::Vector3d> *vn = &grid->v;

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
      }
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	L[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	L[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	L[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
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
  int in;
  vector<Eigen::Vector3d> *vn = &grid->v_update;

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
      }
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	L[ip](0,0) += (*vn)[in][0]*wfd_pn[ip][j][0];
	L[ip](0,1) += (*vn)[in][0]*wfd_pn[ip][j][1];
	L[ip](1,0) += (*vn)[in][1]*wfd_pn[ip][j][0];
	L[ip](1,1) += (*vn)[in][1]*wfd_pn[ip][j][1];
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      L[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
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

void Solid::compute_deformation_gradient()
{
  int in;
  vector<Eigen::Vector3d> *xn = &grid->x;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  Eigen::Vector3d dx;
  Eigen::Matrix3d Ftemp, eye;
  eye.setIdentity();

  if (domain->dimension == 1) {
    for (int ip=0; ip<np_local; ip++){
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*xn)[in] - (*x0n)[in];
	Ftemp(0,0) += dx[0]*wfd_pn[ip][j][0];
      }
      F[ip](0,0) = Ftemp(0,0) +1;
    }
  } else if (domain->dimension == 2) {
    for (int ip=0; ip<np_local; ip++){
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*xn)[in] - (*x0n)[in];
	Ftemp(0,0) += dx[0]*wfd_pn[ip][j][0];
	Ftemp(0,1) += dx[0]*wfd_pn[ip][j][1];
	Ftemp(1,0) += dx[1]*wfd_pn[ip][j][0];
	Ftemp(1,1) += dx[1]*wfd_pn[ip][j][1];
      }
      F[ip] = Ftemp + eye;
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np_local; ip++){
      // F[ip].setZero();
      Ftemp.setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = (*xn)[in] - (*x0n)[in];
	Ftemp(0,0) += dx[0]*wfd_pn[ip][j][0];
	Ftemp(0,1) += dx[0]*wfd_pn[ip][j][1];
	Ftemp(0,2) += dx[0]*wfd_pn[ip][j][2];
	Ftemp(1,0) += dx[1]*wfd_pn[ip][j][0];
	Ftemp(1,1) += dx[1]*wfd_pn[ip][j][1];
	Ftemp(1,2) += dx[1]*wfd_pn[ip][j][2];
	Ftemp(2,0) += dx[2]*wfd_pn[ip][j][0];
	Ftemp(2,1) += dx[2]*wfd_pn[ip][j][1];
	Ftemp(2,2) += dx[2]*wfd_pn[ip][j][2];
      }
      //F[ip].noalias() += eye;
      F[ip] = Ftemp + eye;

    }
  }
}


void Solid::compute_rate_deformation_gradient_TL_APIC()
{
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
  int in;
  vector<Eigen::Vector3d> *x0n = &grid->x0;
  //Eigen::Vector3d *vn = grid->v;
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
  bool status, tl, nh;
  Eigen::Matrix3d U;
  Eigen::Matrix3d eye;
  eye.setIdentity();

  if (update->method_style.compare("tlmpm") == 0) tl = true;
  else tl = false;

  if ((mat->eos==NULL) && (mat->strength==NULL)) nh = true;
  else nh = false;

  for (int ip=0; ip<np_local; ip++){
    
    if (tl) F[ip] += update->dt * Fdot[ip];
    else F[ip] = (eye+update->dt*L[ip]) * F[ip];

    Finv[ip] = F[ip].inverse();
    J[ip] = F[ip].determinant();

    if (J[ip] < 0.0) {
      cout << "Error: J[" << ip << "]<0.0 == " << J[ip] << endl;
      cout << "F[" << ip << "]:" << endl << F[ip] << endl;
      error->all(FLERR,"");
    }

    vol[ip] = J[ip] * vol0[ip];
    rho[ip] = rho0[ip] / J[ip];

    if (!nh) {
      // Only done if not Neo-Hookean:
      if (update->method_style.compare("tlmpm") == 0)
	L[ip] = Fdot[ip] * Finv[ip];
      // else
      //   Fdot[ip] = L[ip]*F[ip];

      status = PolDec(F[ip], R[ip], U, false); // polar decomposition of the deformation gradient, F = R * U
      if (update->method_style.compare("tlmpm") == 0)
	D[ip] = 0.5 * (R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip]);
      else D[ip] = 0.5 * (L[ip] + L[ip].transpose());

      if (!status) {
	cout << "Polar decomposition of deformation gradient failed for particle " << ip << ".\n";
	cout << "F:" << endl << F[ip] << endl;
	cout << "timestep" << endl << update->ntimestep << endl;
	error->all(FLERR,"");
      }
    }

    // strain_increment[ip] = update->dt * D[ip];
  }
}

void Solid::update_stress()
{
  min_inv_p_wave_speed = 1.0e22;
  double pH, plastic_strain_increment;
  Matrix3d eye, sigma_dev, FinvT, PK1;
  bool tl, nh;

  if ((mat->eos!=NULL) && (mat->strength!=NULL)) nh = false;
  else nh = true;

  if (update->method_style.compare("tlmpm") == 0) tl = true;
  else tl = false;

  eye.setIdentity();

  //# pragma omp parallel for
  for (int ip=0; ip<np_local; ip++){
    if (nh) {

      // Neo-Hookean material:
      FinvT = Finv[ip].transpose();
      PK1 = mat->G*(F[ip] - FinvT) + mat->lambda*log(J[ip])*FinvT;
      vol0PK1[ip] = vol0[ip]*PK1;
      sigma[ip] = 1.0/J[ip]*(F[ip]*PK1.transpose());
      strain_el[ip] = 0.5*(F[ip].transpose()*F[ip] - eye);//update->dt * D[ip];

    } else {

      pH = mat->eos->compute_pressure(J[ip], rho[ip], 0, damage[ip]);
      sigma_dev = mat->strength->update_deviatoric_stress(sigma[ip], D[ip], plastic_strain_increment, eff_plastic_strain[ip], eff_plastic_strain_rate[ip], damage[ip]);

      eff_plastic_strain[ip] += plastic_strain_increment;

      // // compute a characteristic time over which to average the plastic strain
      double tav = 1000 * grid->cellsize / mat->signal_velocity;
      eff_plastic_strain_rate[ip] -= eff_plastic_strain_rate[ip] * update->dt / tav;
      eff_plastic_strain_rate[ip] += plastic_strain_increment / tav;
      eff_plastic_strain_rate[ip] = MAX(0.0, eff_plastic_strain_rate[ip]);

      if (mat->damage != NULL)
	mat->damage->compute_damage(damage_init[ip], damage[ip], pH, sigma_dev, eff_plastic_strain_rate[ip], plastic_strain_increment);
      sigma[ip] = -pH*eye + sigma_dev;

      if (damage[ip] > 1e-10) {
	strain_el[ip] = (update->dt*D[ip].trace() + strain_el[ip].trace())/3.0*eye + sigma_dev/(mat->G*(1-damage[ip]));
      } else {
	strain_el[ip] =  (update->dt*D[ip].trace() + strain_el[ip].trace())/3.0*eye;
      }

      if (tl) {
	vol0PK1[ip] = vol0[ip]*J[ip] * (R[ip] * sigma[ip] * R[ip].transpose()) * Finv[ip].transpose();
      }
    }
  }

  double min_h_ratio = 1.0e22;
  double four_third = 1.333333333333333333333333333333333333333;
  for (int ip=0; ip<np_local; ip++){
    min_inv_p_wave_speed = MIN(min_inv_p_wave_speed, rho[ip] / (mat->K + four_third * mat->G));

    min_h_ratio = MIN(min_h_ratio, F[ip](0,0)*F[ip](0,0) + F[ip](0,1)*F[ip](0,1) + F[ip](0,2)*F[ip](0,2));
    min_h_ratio = MIN(min_h_ratio, F[ip](1,0)*F[ip](1,0) + F[ip](1,1)*F[ip](1,1) + F[ip](1,2)*F[ip](1,2));
    min_h_ratio = MIN(min_h_ratio, F[ip](2,0)*F[ip](2,0) + F[ip](2,1)*F[ip](2,1) + F[ip](2,2)*F[ip](2,2));

    if (std::isnan(min_inv_p_wave_speed)) {
      cout << "Error: min_inv_p_wave_speed is nan with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K << ", G=" << mat->G << endl;
      error->all(FLERR, "");
    } else if (min_inv_p_wave_speed < 0.0) {
      cout << "Error: min_inv_p_wave_speed= " << min_inv_p_wave_speed << " with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K << ", G=" << mat->G << endl;
      error->all(FLERR, "");
    }

  }
  min_inv_p_wave_speed = sqrt(min_inv_p_wave_speed);
  dtCFL = MIN(dtCFL, min_inv_p_wave_speed * grid->cellsize * sqrt(min_h_ratio));
  if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "min_inv_p_wave_speed = " << min_inv_p_wave_speed << ", grid->cellsize=" << grid->cellsize << endl;
      error->all(FLERR, "");
  }
}


void Solid::compute_inertia_tensor(string form_function) {

  // int in;
  // Eigen::Vector3d *x0n = grid->x0;
  // Eigen::Vector3d *vn_update = grid->v_update;
  // Eigen::Vector3d dx;

  // if (domain->dimension == 2) {
  //   for (int ip=0; ip<np_local; ip++){
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
  //   for (int ip=0; ip<np_local; ip++){
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

  double cellsizeSqInv = 1.0/(grid->cellsize*grid->cellsize);

  for (int ip=0; ip<np_local; ip++){
    if ( form_function.compare("linear") == 0) {
      // If the form function is linear:
      Di[ip] = 16.0 / 3.0 * cellsizeSqInv * eye;
    } else if ( form_function.compare("quadratic-spline") == 0) {
      // If the form function is a quadratic spline:
      Di[ip] = 4.0 * cellsizeSqInv * eye;
    } else if ( form_function.compare("cubic-spline") == 0) {
      // If the form function is a cubic spline:
      Di[ip] = 3.0 * cellsizeSqInv * eye;
    } else if ( form_function.compare("Bernstein-quadratic") == 0)
      Di[ip] = 12.0 * cellsizeSqInv * eye;
    //cout << "Di[" << ip << "]=\n" << Di[ip] << endl;
  }
}

void Solid::populate(vector<string> args) {

  cout << "Solid delimitated by region ID: " << args[1] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[1]);
  if (iregion == -1) {
    error->all(FLERR, "Error: region ID " + args[1] + " not does not exist.\n");
  }

  if (domain->created==false) {
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
  cout << "proc " << universe->me << "\tsolidsublo=[" << solidsublo[0] << "," << solidsublo[1] << "," << solidsublo[2] << "]\t solidsubhi=["<< solidsubhi[0] << "," << solidsubhi[1] << "," << solidsubhi[2] << "]\n";

  std::vector<double> x2plot, y2plot;
  plt::figure_size(1200, 780);
#endif

  // Calculate total number of particles np_local:
  int nsubx, nsuby, nsubz;
  double delta;
  double hdelta;
  double Lsubx, Lsuby, Lsubz;

  delta = grid->cellsize;
  
  if (grid->nnodes == 0) {
    // The grid will be ajusted to the solid's domain (good for TLMPM).

    // and we need to create the corresponding grid:
    grid->init(solidlo, solidhi);
  }

  Lsubx = solidsubhi[0]-solidsublo[0];
  if (domain->dimension >= 2) Lsuby = solidsubhi[1]-solidsublo[1];
  if (domain->dimension == 3) Lsubz = solidsubhi[2]-solidsublo[2];

  nsubx = (int) (Lsubx/delta);
  while (nsubx*delta <= Lsubx-0.1*delta) nsubx++;
  nsubx++;

  if (domain->dimension >= 2) {
    nsuby = (int) (Lsuby/delta);
    while (nsuby*delta <= Lsuby-0.1*delta) nsuby++;
    nsuby++;
  } else {
    nsuby = 1;
  }

  if (domain->dimension == 3) {
    nsubz = (int) (Lsubz/delta);
    while (nsubz*delta <= Lsubz-0.1*delta) nsubz++;
    nsubz++;
  } else {
    nsubz = 1;
  }

#ifdef DEBUG
  cout << "proc " << universe->me << "\tLsub=[" << Lsubx << "," << Lsuby << "," << Lsubz << "]\t nsub=["<< nsubx << "," << nsuby << "," << nsubz << "]\n";
#endif

  np_local = nsubx*nsuby*nsubz;

  // Create particles:


  cout << "delta = " << delta << endl;

  int l=0;
  double vol_;

  if (domain->dimension == 1) vol_ = delta;
  else if (domain->dimension == 2) vol_ = delta*delta;
  else vol_ = delta*delta*delta;

  double mass_ = mat->rho0 * vol_;

  int np_per_cell = (int) input->parsev(args[2]);

  double *boundlo;
  if (update->method->is_TL)
    boundlo = solidlo;
  else boundlo = domain->boxlo;

  double xi = 0.5;
  double lp = delta;
  int nip = 1;
  vector<double> intpoints;

  if (np_per_cell == 1) {
    // One particle per cell at the center:

    xi = 0.5;
    lp *= 0.5;
    nip = 1;

    intpoints = {0, 0, 0};

  } else if (np_per_cell == 2) {
    // Quadratic elements:

    if (domain->dimension == 1) nip = 2;
    else if (domain->dimension == 2) nip = 4;
    else nip = 8;

    if (nc==0) xi= 0.5/sqrt(3.0);
    else xi = 0.25;

    lp *= 0.25;

    intpoints = {-xi, -xi, -xi,
		 -xi, xi, -xi,
		 xi, -xi, -xi,
		 xi, xi, -xi,
		 -xi, -xi, xi,
		 -xi, xi, xi,
		 xi, -xi, xi,
		 xi, xi, xi};


  } else if (np_per_cell == 3) {
    // Berstein elements:

    if (nc == 0) xi = 0.7746/2;
    else xi = 1.0/3.0;

    lp *= 1.0/6.0;
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
    error->all(FLERR, "Error: solid command 4th argument should be 1 or 2, but " + to_string((int) input->parsev(args[3])) + "received.\n");
  }

  np_local *= nip;
#ifdef DEBUG
  cout << "proc " << universe->me << "\tnp_local=" << np_local << endl;
#endif
  mass_ /= (double) nip;
  vol_ /= (double) nip;

  // Allocate the space in the vectors for np particles:
  grow(np_local);

  int dim = domain->dimension;
  //checkIfInRegion = false;

  // int nsubx0, nsuby0, nsubz0;
  // nsubx0 = (int) (sublo[0] - boundlo[0])/delta;
  // nsuby0 = (int) (sublo[1] - boundlo[1])/delta;
  // nsubz0 = (int) (sublo[2] - boundlo[2])/delta;


  double Loffset[3] = {MAX(0.0, sublo[0] - boundlo[0]),
		       MAX(0.0, sublo[1] - boundlo[1]),
		       MAX(0.0, sublo[2] - boundlo[2])};

  int noffset[3] = {(int) (Loffset[0]/delta),
		    (int) (Loffset[1]/delta),
		    (int) (Loffset[2]/delta)};

  for (int i=0; i<nsubx; i++){
    for (int j=0; j<nsuby; j++){
      for (int k=0; k<nsubz; k++){
	for (int ip=0; ip<nip; ip++) {

	  if (l>=np_local) {
	    cout << "Error in Solid::populate(), exceeding the allocated number of particles.\n";
	    cout << "l = " << l << endl;
	    cout << "np_local = " << np_local << endl;
	    error->all(FLERR, "");
	  }

	  x0[l][0] = x[l][0] = boundlo[0] + delta*(noffset[0] + i + 0.5 + intpoints[3*ip+0]);
	  x0[l][1] = x[l][1] = boundlo[1] + delta*(noffset[1] + j + 0.5 + intpoints[3*ip+1]);
	  if (dim == 3) x0[l][2] = x[l][2] = boundlo[2] + delta*(noffset[2] + k + 0.5 + intpoints[3*ip+2]);
	  else x0[l][2] = x[l][2] = 0;

	  // cout << "x0[" << l << "]=[" << x0[l][0] << "," << x0[l][1] << "," << x0[l][2] << "]\n";;

	  // Check if the particle is inside the region:
	  if (domain->inside_subdomain(x0[l][0], x0[l][1], x0[l][2]) && domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2])==1) {
	    // cout << "Inside\n";

	    if (update->method->is_CPDI && nc != 0) {
	      rp0[dim*l][0] = rp[dim*l][0] = lp;
	      rp0[dim*l][1] = rp[dim*l][1] = 0;
	      rp0[dim*l][2] = rp[dim*l][2] = 0;

	      if (dim >= 2) {
		rp0[dim*l+1][0] = rp[dim*l+1][0] = 0;
		rp0[dim*l+1][1] = rp[dim*l+1][1] = lp;
		rp0[dim*l+1][2] = rp[dim*l+1][2] = 0;

		if (dim == 3) {
		  rp0[dim*l+2][0] = rp[dim*l+1][0] = 0;
		  rp0[dim*l+2][1] = rp[dim*l+1][1] = 0;
		  rp0[dim*l+2][0] = rp[dim*l+1][0] = lp;
		}
	      }
	    }
	    l++;
	  }
	}
      }
    }
  }

  if (np_local > l) {
    grow(l);
  }
  np_local = l; // Adjust np_local to account for the particles outside the domain


  // Determine the total number of particles:
  bigint np_temp = np_local;
  MPI_Allreduce(&np_temp,&np,1,MPI_MPM_BIGINT,MPI_SUM,universe->uworld);
#ifdef DEBUG
  cout << "proc " << universe->me << "\tnp_local=" << np_local << "\tnp=" << np << endl;
#endif

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

#ifdef DEBUG
  cout << "proc " << universe->me << "\tptag0 = " << ptag0 << endl;
#endif


  for (int i=0; i<np_local;i++) {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    mbp[i].setZero();
    v_update[i].setZero();
    vol0[i] = vol[i] = vol_;
    rho0[i] = rho[i] = mat->rho0;
    mass[i] = mass_;
    eff_plastic_strain[i] = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
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
    mask[i] = 1;

    ptag[i] = ptag0 + i + 1;

#ifdef DEBUG
    x2plot.push_back(x0[i][2]);
    y2plot.push_back(x0[i][1]);
#endif
  }

#ifdef DEBUG
  vector<double> xdomain, ydomain;
  vector<double> xsubdomain, ysubdomain;

  xdomain.push_back(domain->boxlo[0]);
  ydomain.push_back(domain->boxlo[1]);

  xdomain.push_back(domain->boxhi[0]);
  ydomain.push_back(domain->boxlo[1]);

  xdomain.push_back(domain->boxhi[0]);
  ydomain.push_back(domain->boxhi[1]);

  xdomain.push_back(domain->boxlo[0]);
  ydomain.push_back(domain->boxhi[1]);

  xdomain.push_back(domain->boxlo[0]);
  ydomain.push_back(domain->boxlo[1]);


  xsubdomain.push_back(sublo[0]);
  ysubdomain.push_back(sublo[1]);

  xsubdomain.push_back(subhi[0]);
  ysubdomain.push_back(sublo[1]);

  xsubdomain.push_back(subhi[0]);
  ysubdomain.push_back(subhi[1]);

  xsubdomain.push_back(sublo[0]);
  ysubdomain.push_back(subhi[1]);

  xsubdomain.push_back(sublo[0]);
  ysubdomain.push_back(sublo[1]);

  plt::plot(xdomain, ydomain, "b-");
  plt::plot(xsubdomain, ysubdomain, "r-");
  plt::plot(x2plot, y2plot, ".");
  plt::axis("equal");
  plt::save("debug-proc_" + to_string(universe->me) + ".png");
  plt::close();
#endif
}

void Solid::update_particle_domain() {
  int dim = domain->dimension;

  for (int ip=0; ip<np_local; ip++){
    rp[dim*ip] = F[ip]*rp0[dim*ip];
    if (dim >= 2) rp[dim*ip+1] = F[ip]*rp0[dim*ip+1];
    if (dim == 3) rp[dim*ip+2] = F[ip]*rp0[dim*ip+2];
  }
}
