#include "mpm.h"
#include "solid.h"
#include "material.h"
#include "memory.h"
#include <vector>

using namespace std;


Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new solid with ID: " << args[0] << endl;
  id = args[0];

  np = 0;

  x = x0 = NULL;

  vol = vol0 = NULL;
  mass = NULL;
  mask = NULL;

  eos = NULL;
  grid = new Grid(mpm);
}

Solid::~Solid()
{
  memory->destroy(x0);
  memory->destroy(x);
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
  x0 = memory->grow(x0, np, 3, str);

  str = "solid-" + id + ":x";
  cout << "Growing " << str << endl;
  x = memory->grow(x, np, 3, str);

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
