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

  int np, nx, ny, nz;
  nx = input->parse(args[3]);
  ny = input->parse(args[4]);
  nz = input->parse(args[5]);
  np = nx*ny*nz;

  grow(np);

  vector<double> limits = domain->regions[iregion]->limits();
  double delta_x, delta_y, delta_z;
  delta_x = (limits[1]-limits[0])/((float) nx);
  delta_y = (limits[3]-limits[2])/((float) ny);
  delta_z = (limits[5]-limits[4])/((float) nz);

  cout << "deltas = " << delta_x << "\t"<< delta_y << "\t"<< delta_z << "\t" << endl;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	
      }
    }
  }
}


SolBlock::~SolBlock()
{

}
