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

  solidlo[0] = limits[0];
  solidhi[0] = limits[1];
  solidlo[1] = limits[2];
  solidhi[1] = limits[3];
  solidlo[2] = limits[4];
  solidhi[2] = limits[5];

  double delta_x, delta_y, delta_z;
  double hdelta_x, hdelta_y, hdelta_z;

  delta_x = (solidhi[0]-solidlo[0])/((float) nx);
  delta_y = (solidhi[1]-solidlo[1])/((float) ny);
  delta_z = (solidhi[2]-solidlo[2])/((float) nz);

  hdelta_x = 0.5*delta_x;
  hdelta_y = 0.5*delta_y;
  hdelta_z = 0.5*delta_z;

  cout << "deltas = " << delta_x << "\t"<< delta_y << "\t"<< delta_z << "\t" << endl;

  int l=0;
  double vol_ = delta_x*delta_y*delta_z;
  double mass_ = eos->rho0() * vol_;

  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	x0[l][0] = x[l][0] = solidlo[0] + delta_x*(i+0.5);
	x0[l][1] = x[l][1] = solidlo[1] + delta_y*(j+0.5);
	x0[l][2] = x[l][2] = solidlo[2] + delta_z*(k+0.5);


	l++;
      }
    }
  }

  for (int i=0; i<np;i++) {
    v[i].setZero();
    f[i].setZero();
    b[i].setZero();
    v_update[i].setZero();
    vol0[i] = vol[i] = vol_;
    rho0[i] = rho[i] = eos->rho0();
    mass[i] = mass_;
    sigma[i].setZero();
    PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setZero();
    U[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();
    strain_increment[i].setZero();

    J[i] = 1;
  }

  if (l!=np) {
    cout << "Error l=" << l << " != np=" << np << endl;
    exit(1);
  }
}


SolBlock::~SolBlock()
{

}
