#include "mpm.h"
#include "solid.h"
#include "material.h"
#include "memory.h"
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new solid with ID: " << args[0] << endl;
  id = args[0];

  np = 0;

  x = x0 = NULL;
  v = v_update = NULL;

  a = NULL;

  sigma = PK1 = L = F = R = U = Fdot = strain_increment = NULL;

  b = f = NULL;

  J = NULL;

  vol = vol0 = NULL;
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
  if (Fdot!=NULL) delete Fdot;
  if (strain_increment!=NULL) delete strain_increment;

  memory->destroy(J);
  memory->destroy(vol);
  memory->destroy(vol0);
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
