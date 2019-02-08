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

  b = f = NULL;

  vol = vol0 = NULL;
  mass = NULL;
  mask = NULL;

  eos = NULL;
  grid = new Grid(mpm);
}

Solid::~Solid()
{
  delete [] x0;
  delete [] x;
  delete [] v;
  delete [] v_update;
  delete [] b;
  delete [] f;

  memory->destroy(vol);
  memory->destroy(vol0);
  memory->destroy(mass);
  memory->destroy(mask);

  delete eos;
  delete grid;
}


void Solid::init()
{
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

    grid->init(*it); // set the grid cellsize

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
  if (x0 == NULL)
    x0 = new Eigen::Vector3d[np];
  else {
    cout << "Error: x0 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":x";
  cout << "Growing " << str << endl;
  if (x == NULL)
    x = new Eigen::Vector3d[np];
  else {
    cout << "Error: x already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":v";
  cout << "Growing " << str << endl;
  if (v == NULL)
    v = new Eigen::Vector3d[np];
  else {
    cout << "Error: v already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":v_update";
  cout << "Growing " << str << endl;
  if (v_update == NULL)
    v_update = new Eigen::Vector3d[np];
  else {
    cout << "Error: v_update already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":b";
  cout << "Growing " << str << endl;
  if (b == NULL)
    b = new Eigen::Vector3d[np];
  else {
    cout << "Error: b already exists, I don't know how to grow it!\n";
    exit(1);
  }

  str = "solid-" + id + ":f";
  cout << "Growing " << str << endl;
  if (f == NULL)
    f = new Eigen::Vector3d[np];
  else {
    cout << "Error: f already exists, I don't know how to grow it!\n";
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
}
