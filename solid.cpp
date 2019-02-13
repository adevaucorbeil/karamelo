#include "mpm.h"
#include "solid.h"
#include "material.h"
#include "memory.h"
#include "update.h"
#include "domain.h"
#include <vector>
#include <Eigen/Eigen>
#include "mpm_math.h"

using namespace std;
using namespace Eigen;
using namespace MPM_Math;

Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new solid with ID: " << args[0] << endl;
  id = args[0];

  np = 0;

  x = x0 = NULL;
  v = v_update = NULL;

  a = NULL;

  sigma = PK1 = L = F = R = U = Finv = Fdot = strain_increment = NULL;

  b = f = NULL;

  J = NULL;

  vol = vol0 = NULL;
  rho = rho0 = NULL;
  mass = NULL;
  mask = NULL;

  eos = NULL;
  grid = new Grid(mpm);

  numneigh_pn = numneigh_np = NULL;

  neigh_pn = neigh_np = NULL;

  wf_pn = wf_np = NULL;
  wfd_pn = wfd_np = NULL;
}

Solid::~Solid()
{
  if (x0!=NULL) delete x0;
  if (x!=NULL) delete x;
  if (v!=NULL) delete v;
  if (v_update!=NULL) delete v_update;
  if (a!=NULL) delete a;
  if (b!=NULL) delete b;
  if (f!=NULL) delete f;
  if (sigma!=NULL) delete sigma;
  if (PK1!=NULL) delete PK1;
  if (L!=NULL) delete L;
  if (F!=NULL) delete F;
  if (R!=NULL) delete R;
  if (U!=NULL) delete U;
  if (Finv!=NULL) delete Finv;
  if (Fdot!=NULL) delete Fdot;
  if (strain_increment!=NULL) delete strain_increment;

  memory->destroy(J);
  memory->destroy(vol);
  memory->destroy(vol0);
  memory->destroy(rho);
  memory->destroy(rho0);
  memory->destroy(mass);
  memory->destroy(mask);

  delete eos;
  delete grid;

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

  grid->init(solidlo, solidhi);

  if (np == 0) {
    cout << "Error: solid does not have any particles" << endl;
  } else {
      bigint nnodes = grid->nnodes;

      numneigh_pn = new int[np]();
      neigh_pn = new vector<int>[np];
      wf_pn = new vector<double>[np];
      wfd_pn = new vector< array<double,3> >[np];

      if (nnodes) {
	numneigh_np = new int[nnodes]();
	neigh_np = new vector<int>[nnodes];
	wf_np = new vector<double>[nnodes];
	wfd_np = new vector< array<double,3> >[nnodes];
      }
  }
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
  if (args->end() < it+2) {
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it) {
    int iEOS = material->find_EOS(*it);

    if (iEOS == -1){
      cout << "Error: could not find EOS named " << *it << endl;
      exit(1);
    }

    eos = material->EOSs[iEOS]; // point eos to the right EOS class

    it++;

    grid->setup(*it); // set the grid cellsize

    it++;

    if (it != args->end()) {
      cout << "Error: too many arguments" << endl;
      exit(1);
    }
  }
}


void Solid::grow(int nparticles){
  np = nparticles;

  string str;
  str = "solid-" + id + ":x0";
  cout << "Growing " << str << endl;
  if (x0 == NULL) x0 = new Eigen::Vector3d[np];
  else {
    cout << "Error: x0 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":x";
  cout << "Growing " << str << endl;
  if (x == NULL) x = new Eigen::Vector3d[np];
  else {
    cout << "Error: x already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":v";
  cout << "Growing " << str << endl;
  if (v == NULL) v = new Eigen::Vector3d[np];
  else {
    cout << "Error: v already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":v_update";
  cout << "Growing " << str << endl;
  if (v_update == NULL) v_update = new Eigen::Vector3d[np];
  else {
    cout << "Error: v_update already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":a";
  cout << "Growing " << str << endl;
  if (a == NULL) a = new Eigen::Vector3d[np];
  else {
    cout << "Error: a already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":b";
  cout << "Growing " << str << endl;
  if (b == NULL) b = new Eigen::Vector3d[np];
  else {
    cout << "Error: b already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":f";
  cout << "Growing " << str << endl;
  if (f == NULL) f = new Eigen::Vector3d[np];
  else {
    cout << "Error: f already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (sigma == NULL) sigma = new Eigen::Matrix3d[np];
  else {
    cout << "Error: sigma already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (PK1 == NULL) PK1 = new Eigen::Matrix3d[np];
  else {
    cout << "Error: PK1 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (L == NULL) L = new Eigen::Matrix3d[np];
  else {
    cout << "Error: L already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (F == NULL) F = new Eigen::Matrix3d[np];
  else {
    cout << "Error: F already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (R == NULL) R = new Eigen::Matrix3d[np];
  else {
    cout << "Error: R already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (U == NULL) U = new Eigen::Matrix3d[np];
  else {
    cout << "Error: U already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (Finv == NULL) Finv = new Eigen::Matrix3d[np];
  else {
    cout << "Error: Finv already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (Fdot == NULL) Fdot = new Eigen::Matrix3d[np];
  else {
    cout << "Error: Fdot already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (strain_increment == NULL) strain_increment = new Eigen::Matrix3d[np];
  else {
    cout << "Error: strain_increment already exists, I don't know how to grow it!\n";
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

  str = "solid-" + id + ":mask";
  cout << "Growing " << str << endl;
  mask = memory->grow(mask, np, str);

  for (int i=0; i<np; i++) mask[i] = 1;

  str = "solid-" + id + ":J";
  cout << "Growing " << str << endl;
  J = memory->grow(J, np, str);
}

void Solid::compute_mass_nodes()
{
  int ip;
  for (int in=0; in<grid->nnodes; in++){
    grid->mass[in] = 0;
    for (int j=0; j<numneigh_np[in];j++){
      ip = neigh_np[in][j];
      grid->mass[in] += wf_np[in][j] * mass[ip];
    }
  }
  return;
}

void Solid::compute_velocity_nodes()
{
  Eigen::Vector3d *vn = grid->v;
  double *massn = grid->mass;
  int ip;
  
  for (int in=0; in<grid->nnodes; in++) {
    vn[in].setZero();
    if (massn[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	vn[in] += (wf_np[in][j] * mass[ip] / massn[in]) * v[ip];
      }
    }
  }
}

void Solid::compute_external_forces_nodes()
{
  Eigen::Vector3d *bn = grid->b;
  double *massn = grid->mass;
  int ip;
  
  for (int in=0; in<grid->nnodes; in++) {
    bn[in].fill(0);
    if (massn[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	bn[in] += (wf_np[in][j] * mass[ip] / massn[in]) * b[ip];
      }
    }
  }
}

void Solid::compute_internal_forces_nodes()
{
  Eigen::Vector3d *fn = grid->f;
  double *massn = grid->mass;
  int ip;
  
  for (int in=0; in<grid->nnodes; in++) {
    fn[in].fill(0);
    if (massn[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	fn[in] += (wf_np[in][j] * mass[ip] / massn[in]) * f[ip];
      }
    }
  }
}

void Solid::compute_particle_velocities()
{
  Eigen::Vector3d *vn_update = grid->v_update;
  int in;

  for (int ip=0; ip<np; ip++){
    v_update[ip].fill(0);
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * vn_update[in];
    }
  }
}

void Solid::compute_particle_acceleration()
{
  double inv_dt = 1.0/update->dt;
  
  Eigen::Vector3d *vn_update = grid->v_update;
  Eigen::Vector3d *vn = grid->v;

  int in;

  for (int ip=0; ip<np; ip++){
    a[ip].fill(0);
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      a[ip] +=inv_dt * wf_pn[ip][j] * (vn_update[in] - vn[j]);
    }
    f[ip] = a[ip] / mass[ip];
  }
}

void Solid::update_particle_position()
{
  for (int ip=0; ip<np; ip++) {
    x[ip] += update->dt*v_update[ip];
  }
}

void Solid::update_particle_velocities(double FLIP)
{
  for (int ip=0; ip<np; ip++) {
    v[ip] = (1 - FLIP) * v_update[ip] + FLIP*(v[ip] + update->dt*a[ip]);
  }
}

void Solid::compute_rate_deformation_gradient()
{
  int in;
  Eigen::Vector3d *vn = grid->v;

  if (domain->dimension == 2) {
    for (int ip=0; ip<np; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	Fdot[ip](0,0) += vn[in][0]*wfd_pn[ip][j][0];
	Fdot[ip](0,1) += vn[in][0]*wfd_pn[ip][j][1];
	Fdot[ip](1,0) += vn[in][1]*wfd_pn[ip][j][0];
	Fdot[ip](1,1) += vn[in][1]*wfd_pn[ip][j][1];
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	Fdot[ip](0,0) += vn[in][0]*wfd_pn[ip][j][0];
	Fdot[ip](0,1) += vn[in][0]*wfd_pn[ip][j][1];
	Fdot[ip](0,2) += vn[in][0]*wfd_pn[ip][j][2];
	Fdot[ip](1,0) += vn[in][1]*wfd_pn[ip][j][0];
	Fdot[ip](1,1) += vn[in][1]*wfd_pn[ip][j][1];
	Fdot[ip](1,2) += vn[in][1]*wfd_pn[ip][j][2];
	Fdot[ip](2,0) += vn[in][2]*wfd_pn[ip][j][0];
	Fdot[ip](2,1) += vn[in][2]*wfd_pn[ip][j][1];
	Fdot[ip](2,2) += vn[in][2]*wfd_pn[ip][j][2];
      }
    }
  }
}

void Solid::update_deformation_gradient()
{
  bool status;
  Eigen::Matrix3d U;

  for (int ip=0; ip<np; ip++){
    F[ip] += update->dt * Fdot[ip];
    J[ip] = F[ip].determinant();
    vol[ip] = J[ip] * vol0[ip];
    rho[ip] = rho0[ip] / J[ip];
    Finv[ip] = F[ip].inverse();
    L[ip] = Fdot[ip] * Finv[ip];

    status = PolDec(F[ip], R[ip], U, false); // polar decomposition of the deformation gradient, F = R * U

    if (!status) {
      cout << "Polar decomposition of deformation gradient failed for particle " << ip << ".\n";
      exit(1);
    }

    strain_increment[ip] = 0.5*update->dt * R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip];
  }
}

void Solid::update_stress()
{
  for (int ip=0; ip<np; ip++){
    eos->update_stress(sigma[ip], strain_increment[ip], J[ip]);
    PK1[ip] = J[ip] * (R[ip] * sigma[ip] * R[ip].transpose()) * Finv[ip].transpose();
  }

}
