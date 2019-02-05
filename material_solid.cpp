#include <iostream>
#include "material_solid.h"
#include "input.h"

using namespace std;


MatSol::MatSol(MPM *mpm, vector<string> args) : Mat(mpm, args)
{
  cout << "Initiate MatSol" << endl;

  if (args.size()<8) {
    cout << "Error: material command not enough arguments" << endl;
    exit(1);
  }
  options(&args, args.begin()+8);

  rho0 = 0;
  cp = 0;
  K0 = 0;
  G0 = 0;
  thermal_conductivity = 0;
  thermal_diffusivity = 0;

  
}


MatSol::~MatSol()
{

}

