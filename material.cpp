#include "mpm.h"
#include "material.h"
#include <vector>

using namespace std;


Material::Material(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new material with ID: " << args[0] << endl;
  id = args[0];

  rho0 = 0;
  cp = 0;
  K0 = 0;
  G0 = 0;
  thermal_conductivity = 0;
  thermal_diffusivity = 0;
}

Material::~Material()
{
}


void Material::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In material::options()" << endl;
  if (args->end() < it) {
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

