#include "mpm.h"
#include "solid.h"
#include "memory.h"
#include "update.h"
#include "domain.h"
#include <vector>
#include <Eigen/Eigen>
#include "mpm_math.h"
#include <math.h>

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

  sigma = PK1 = L = F = R = U = D = Finv = Fdot = Di = NULL;

  b = f = NULL;

  J = NULL;

  vol = vol0 = NULL;
  rho = rho0 = NULL;
  mass = NULL;
  eff_plastic_strain = NULL;
  eff_plastic_strain_rate = NULL;
  damage = NULL;
  damage_init = NULL;
  mask = NULL;

  mat = NULL;
  grid = new Grid(mpm);

  numneigh_pn = numneigh_np = NULL;

  neigh_pn = neigh_np = NULL;

  wf_pn = wf_np = NULL;
  wfd_pn = wfd_np = NULL;

  dtCFL = 1.0e22;
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
  if (D!=NULL) delete D;
  if (Finv!=NULL) delete Finv;
  if (Fdot!=NULL) delete Fdot;
  if (Di!=NULL) delete Di;

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
  memory->destroy(mask);

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
      wfd_pn = new vector< Vector3d >[np];

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
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it) {
    int iMat = material->find_material(*it);

    if (iMat == -1){
      cout << "Error: could not find material named " << *it << endl;
      exit(1);
    }

    mat = &material->materials[iMat]; // point mat to the right material

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

  if (D == NULL) D = new Eigen::Matrix3d[np];
  else {
    cout << "Error: D already exists, I don't know how to grow it!\n";
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

  if (Di == NULL) Di = new Eigen::Matrix3d[np];
  else {
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
	// if (update->ntimestep>450)
	//   if (in==0)
	//     cout << "compute_velocity_nodes: in=" << in << ", ip=" << ip << ", j= " << j << ", wf_np=" << wf_np[in][j] << ", vn=[" << vn[in][0] << "," << vn[in][1] << "," << vn[in][2] << "], vp=[" << v[ip][0] << "," << v[ip][1] << "," << v[ip][2] << "]" << endl;
      }
    }
  }
}

void Solid::compute_velocity_nodes_APIC()
{
  Eigen::Vector3d *x0n = grid->x0;
  Eigen::Vector3d *vn = grid->v;
  double *massn = grid->mass;
  int ip;
  
  for (int in=0; in<grid->nnodes; in++) {
    vn[in].setZero();
    if (massn[in] > 0) {
      for (int j=0; j<numneigh_np[in];j++){
	ip = neigh_np[in][j];
	vn[in] += (wf_np[in][j] * mass[ip] / massn[in]) * (v[ip] + Fdot[ip]*(x0n[in] - x0[ip]));
	// if (in==32)
	//   cout << "compute_velocity_nodes: in=" << in << ", ip=" << ip << ", j= " << j << ", wf_np=" << wf_np[in][j] << ", vn=[" << vn[in][0] << "," << vn[in][1] << "," << vn[in][2] << "], vp=[" << v[ip][0] << "," << v[ip][1] << "," << v[ip][2] << "]" << endl;
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
    bn[in].setZero();
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
    fn[in].setZero();
    for (int j=0; j<numneigh_np[in];j++){
      ip = neigh_np[in][j];
      fn[in] -= vol0[ip] * (PK1[ip] * wfd_np[in][j]);
    }
  }
}

void Solid::compute_particle_velocities()
{
  Eigen::Vector3d *vn_update = grid->v_update;
  int in;

  for (int ip=0; ip<np; ip++){
    v_update[ip].setZero();
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      v_update[ip] += wf_pn[ip][j] * vn_update[in];
	// if (update->ntimestep>450)
	//   if (ip==0)
	//     cout << "compute_particle_velocities: ip=" << ip << ", in=" << in << ", j= " << j << ", wf_pn=" << wf_pn[ip][j] << ", vp_update=[" << v_update[ip][0] << "," << v_update[ip][1] << "," << v_update[ip][2] << "], vn_update=[" << vn_update[in][0] << "," << vn_update[in][1] << "," << vn_update[in][2] << "]" << endl;
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
    a[ip].setZero();
    for (int j=0; j<numneigh_pn[ip]; j++){
      in = neigh_pn[ip][j];
      a[ip] += inv_dt * wf_pn[ip][j] * (vn_update[in] - vn[in]);
	// if (update->ntimestep>450)
	//   if (ip==0)
	//     cout << "compute_particle_acceleration: ip=" << ip  << ", in=" << in << ", ap = [" << a[ip][0]  << "," <<  a[ip][1] << "," <<  a[ip][2] << "], vn_update = [" <<  vn_update[in][0] << "," << vn_update[in][1]  << "," << vn_update[in][2] << "], vn = [" << vn[in][0]  << "," << vn[in][1]  << "," <<  vn[in][2] << "], wf=" <<  wf_pn[ip][j] << ", inv_dt=" << inv_dt << endl;
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
	// if (update->ntimestep>450)
	//   if (ip==0)
	//     cout << "update_particle_velocities: ip=" << ip << ", FLIP=" << FLIP << ", v[ip]=[" << v[ip][0] << "," << v[ip][1] << "," << v[ip][2] << "], dt=" << update->dt << ", ap=[" << a[ip][0] << "," << a[ip][1] << "," << a[ip][2] << "]" << endl;
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
	// if (update->ntimestep > 450)
	//   if (ip==0) cout << "compute_rate_deformation_gradient: ip=" << ip << ", in=" << in << ", vn=[" << vn[in][0] << "," << vn[in][1] << "," << vn[in][2] << "], wfd_pn=[" << wfd_pn[ip][j][0] << "," << wfd_pn[ip][j][1] << "," << wfd_pn[ip][j][2] << "], Fdot=[[" << Fdot[ip](0,0) << "," << Fdot[ip](0,1) << "," << Fdot[ip](0,2) << "],[" << Fdot[ip](1,0) << "," << Fdot[ip](1,1) << "," << Fdot[ip](1,2) << "],[" << Fdot[ip](2,0) << "," << Fdot[ip](2,1) << "," << Fdot[ip](2,2) << "]]" << endl;
      }
    }
  }
}
void Solid::compute_deformation_gradient()
{
  int in;
  Eigen::Vector3d *xn = grid->x;
  Eigen::Vector3d *x0n = grid->x0;
  Eigen::Vector3d dx;
  Eigen::Matrix3d eye;
  eye.setIdentity();

  if (domain->dimension == 2) {
    for (int ip=0; ip<np; ip++){
      F[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = xn[in] - x0n[in];
	F[ip](0,0) += dx[0]*wfd_pn[ip][j][0];
	F[ip](0,1) += dx[0]*wfd_pn[ip][j][1];
	F[ip](1,0) += dx[1]*wfd_pn[ip][j][0];
	F[ip](1,1) += dx[1]*wfd_pn[ip][j][1];
      }
      F[ip].noalias() += eye;
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np; ip++){
      F[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = xn[in] - x0n[in];
	F[ip](0,0) += dx[0]*wfd_pn[ip][j][0];
	F[ip](0,1) += dx[0]*wfd_pn[ip][j][1];
	F[ip](0,2) += dx[0]*wfd_pn[ip][j][2];
	F[ip](1,0) += dx[1]*wfd_pn[ip][j][0];
	F[ip](1,1) += dx[1]*wfd_pn[ip][j][1];
	F[ip](1,2) += dx[1]*wfd_pn[ip][j][2];
	F[ip](2,0) += dx[2]*wfd_pn[ip][j][0];
	F[ip](2,1) += dx[2]*wfd_pn[ip][j][1];
	F[ip](2,2) += dx[2]*wfd_pn[ip][j][2];
      }
      F[ip].noalias() += eye;

    }
  }
}


void Solid::compute_rate_deformation_gradient_APIC()
{
  int in;
  Eigen::Vector3d *x0n = grid->x0;
  //Eigen::Vector3d *vn = grid->v;
  Eigen::Vector3d *vn = grid->v_update;
  Eigen::Vector3d dx;

  if (domain->dimension == 2) {
    for (int ip=0; ip<np; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = x0n[in] - x0[ip];
	Fdot[ip](0,0) += vn[in][0]*dx[0]*wf_pn[ip][j];
	Fdot[ip](0,1) += vn[in][0]*dx[1]*wf_pn[ip][j];
	Fdot[ip](1,0) += vn[in][1]*dx[0]*wf_pn[ip][j];
	Fdot[ip](1,1) += vn[in][1]*dx[1]*wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
    }
  } else if (domain->dimension == 3) {
    for (int ip=0; ip<np; ip++){
      Fdot[ip].setZero();
      for (int j=0; j<numneigh_pn[ip]; j++){
	in = neigh_pn[ip][j];
	dx = x0n[in] - x0[ip];
	Fdot[ip](0,0) += vn[in][0]*dx[0]*wf_pn[ip][j];
	Fdot[ip](0,1) += vn[in][0]*dx[1]*wf_pn[ip][j];
	Fdot[ip](0,2) += vn[in][0]*dx[2]*wf_pn[ip][j];
	Fdot[ip](1,0) += vn[in][1]*dx[0]*wf_pn[ip][j];
	Fdot[ip](1,1) += vn[in][1]*dx[1]*wf_pn[ip][j];
	Fdot[ip](1,2) += vn[in][1]*dx[2]*wf_pn[ip][j];
	Fdot[ip](2,0) += vn[in][2]*dx[0]*wf_pn[ip][j];
	Fdot[ip](2,1) += vn[in][2]*dx[1]*wf_pn[ip][j];
	Fdot[ip](2,2) += vn[in][2]*dx[2]*wf_pn[ip][j];
      }
      Fdot[ip] *= Di[ip];
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
    if (J[ip] < 0.0) {
      cout << "Error: J[" << ip << "]<0.0 == " << J[ip] << endl;
      cout << "F[" << ip << "]:" << endl << F[ip] << endl;
	exit(1);
    }
    vol[ip] = J[ip] * vol0[ip];
    rho[ip] = rho0[ip] / J[ip];
    Finv[ip] = F[ip].inverse();
    L[ip] = Fdot[ip] * Finv[ip];

    status = PolDec(F[ip], R[ip], U, false); // polar decomposition of the deformation gradient, F = R * U
    D[ip] = 0.5 * (R[ip].transpose() * (L[ip] + L[ip].transpose()) * R[ip]);

    if (!status) {
      cout << "Polar decomposition of deformation gradient failed for particle " << ip << ".\n";
      cout << "F:" << endl << F[ip] << endl;
      cout << "timestep" << endl << update->ntimestep << endl;
      exit(1);
    }

    // strain_increment[ip] = update->dt * D[ip];
  }
}

void Solid::update_stress()
{
  min_inv_p_wave_speed = 1.0e22;
  double pH, plastic_strain_increment;
  Matrix3d eye, sigma_dev, FinvT;

  eye.setIdentity();

  for (int ip=0; ip<np; ip++){
    if ((mat->eos!=NULL) && (mat->strength!=NULL)) {
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

      PK1[ip] = J[ip] * (R[ip] * sigma[ip] * R[ip].transpose()) * Finv[ip].transpose();

    } else {
      // Neo-Hookean material:
      FinvT = Finv[ip].transpose();
      PK1[ip] = mat->G*(F[ip] - FinvT) + mat->lambda*log(J[ip])*FinvT;
      sigma[ip] = 1.0/J[ip]*(PK1[ip]*F[ip].transpose());
    }

    
    min_inv_p_wave_speed = MIN(min_inv_p_wave_speed, rho[ip] / (mat->K + 4.0/3.0 * mat->G));
    if (std::isnan(min_inv_p_wave_speed)) {
      cout << "Error: min_inv_p_wave_speed is nan with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K << ", G=" << mat->G << endl;
      exit(1);
    } else if (min_inv_p_wave_speed < 0.0) {
      cout << "Error: min_inv_p_wave_speed= " << min_inv_p_wave_speed << " with ip=" << ip << ", rho[ip]=" << rho[ip] << ", K=" << mat->K << ", G=" << mat->G << endl;
      exit(1);
    }

  }
  min_inv_p_wave_speed = sqrt(min_inv_p_wave_speed);
  dtCFL = MIN(dtCFL, min_inv_p_wave_speed * grid->cellsize);
  if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "min_inv_p_wave_speed = " << min_inv_p_wave_speed << ", grid->cellsize=" << grid->cellsize << endl;
      exit(1);
  }
}

void Solid::compute_inertia_tensor(string form_function) {

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

  Eigen::Matrix3d eye;
  eye.setIdentity();

  double cellsizeSqInv = 1.0/(grid->cellsize*grid->cellsize);

  for (int ip=0; ip<np; ip++){
    if ( form_function.compare("linear") == 0) {
      // If the form function is linear:
      Di[ip] = 16.0 / 3.0 * cellsizeSqInv * eye;
    } else if ( form_function.compare("quadratic-spline") == 0) {
      // If the form function is a quadratic spline:
      Di[ip] = 4.0 * cellsizeSqInv * eye;
    } else if ( form_function.compare("cubic-spline") == 0) {
      // If the form function is a cubic spline:
      Di[ip] = 3.0 * cellsizeSqInv * eye;
    }
    //cout << "Di[" << ip << "]=\n" << Di[ip] << endl;
  }
}
