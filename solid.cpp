#include "mpm.h"
#include "solid.h"
#include "material.h"
#include <vector>

using namespace std;


Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new solid with ID: " << args[0] << endl;
  id = args[0];

  eos = NULL;
  grid = new Grid(mpm);
}

Solid::~Solid()
{
  delete eos;
  delete grid;
}


void Solid::init()
{
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
  if (args->end() < it+1) {
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

    if (it != args->end()) {
      cout << "Error: too many arguments" << endl;
      exit(1);
    }
  }
}


void Solid::grow(int np){
  nparticles = np;
  x0.reserve(nparticles);
  x.reserve(nparticles);
  vol0.reserve(nparticles);
  vol.reserve(nparticles);
  mass.reserve(nparticles);
}
