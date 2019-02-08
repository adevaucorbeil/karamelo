#include <iostream>
#include "domain.h"
#include "solid_block.h"
#include "input.h"

using namespace std;


SolBlock::SolBlock(MPM *mpm, vector<string> args) : Solid(mpm, args)
{
  cout << "Initiate SolBlock" << endl;

  if (args.size() < domain->dimension + 3) {
    cout << "Error: solid command not enough arguments" << endl;
    exit(1);
  }
  options(&args, args.begin()+domain->dimension + 3);

  cout << "Solid delimitated by region ID: " << args[2] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1) {
    cout << "Error: region ID " << args[2] << " not does not exist" << endl;
    exit(1);
  }

  // Calculate total number of particles np:
  int np, nx, ny, nz;
  nx = input->parse(args[3]);
  ny = input->parse(args[4]);
  nz = input->parse(args[5]);
  np = nx*ny*nz;

  // Allocate the space in the vectors for np particles:
  grow(np);

  // Create particles:
  vector<double> limits = domain->regions[iregion]->limits();
  double delta_x, delta_y, delta_z;
  delta_x = (limits[1]-limits[0])/((float) nx);
  delta_y = (limits[3]-limits[2])/((float) ny);
  delta_z = (limits[5]-limits[4])/((float) nz);

  cout << "deltas = " << delta_x << "\t"<< delta_y << "\t"<< delta_z << "\t" << endl;

  int l=0;
  double vol_ = delta_x*delta_y*delta_z;
  double mass_ = eos->rho0() * vol_;

  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	x0[l][0] = x[l][0] = limits[0] + delta_x*(i+0.5);
	x0[l][1] = x[l][1] = limits[2] + delta_y*(j+0.5);
	x0[l][2] = x[l][2] = limits[4] + delta_z*(k+0.5);
	vol0[l] = vol[l] = vol_;
	mass[l] = mass_;
	l++;
      }
    }
  }
}


SolBlock::~SolBlock()
{

}
