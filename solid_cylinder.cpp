#include <iostream>
#include "domain.h"
#include "solid_cylinder.h"
#include "input.h"
#include "var.h"
#include "grid.h"

using namespace std;

SolCylinder::SolCylinder(MPM *mpm, vector<string> args) : Solid(mpm, args)
{
  cout << "Initiate SolCylinder" << endl;

  if (args.size() < 4) {
    cout << "Error: solid command not enough arguments" << endl;
    exit(1);
  }
  options(&args, args.begin()+4);

  cout << "Solid delimitated by region ID: " << args[2] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1) {
    cout << "Error: region ID " << args[2] << " not does not exist" << endl;
    exit(1);
  }

  vector<double> limits = domain->regions[iregion]->limits();

  solidlo[0] = limits[0];
  solidhi[0] = limits[1];
  solidlo[1] = limits[2];
  solidhi[1] = limits[3];
  solidlo[2] = limits[4];
  solidhi[2] = limits[5];

  if (grid->nnodes == 0) grid->init(solidlo, solidhi);

  // Calculate total number of particles np:
  int nx, ny, nz;
  double delta;
  double hdelta;

  delta = grid->cellsize;///((int) input->parsev(args[3]));

  double Lx = solidhi[0]-solidlo[0];
  double Ly = solidhi[1]-solidlo[1];

  nx = (int) (Lx/delta);
  while (nx * delta <= Lx - delta) {
    nx++;
  }

  ny = (int) (Ly/delta);
  while (ny * delta <= Ly - delta) {
    ny++;
  }

  if (domain->dimension == 3) {
    double Lz = solidhi[2]-solidlo[2];
    nz = (int) (Lz/delta);

    while (nz * delta <= Lz - delta) {
      nz++;
    }
   } else {
    nz = 1;
   }

  np = nx*ny*nz;


  // Create particles:


  cout << "delta = " << delta << endl;

  int l=0;
  double vol_;
  if (domain->dimension == 3) vol_ = delta*delta*delta;
  else vol_ = delta*delta;

  double mass_ = mat->rho0 * vol_;

  if ((int) input->parsev(args[3]) == 1) {
    // One particle per cell at the center:

    // Allocate the space in the vectors for np particles:
    grow(np);

    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
	for (int k=0; k<nz; k++){
	  x0[l][0] = x[l][0] = solidlo[0] + delta*(i+0.5);
	  x0[l][1] = x[l][1] = solidlo[1] + delta*(j+0.5);
	  if (domain->dimension == 3) x0[l][2] = x[l][2] = solidlo[2] + delta*(k+0.5);
	  else x0[l][2] = x[l][2] = 0;
	  // Check if the particle is inside cylinder:
	  if (domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2])==1) l++;
	}
      }
    }
  } else if ((int) input->parsev(args[3]) == 2) {
    // Quadratic elements:

    np *= 8;
    mass_ /= 8.0;
    vol_ /= 8.0;
    // Allocate the space in the vectors for np particles:
    grow(np);

    double half_Sqrt_three_inv = 0.5/sqrt(3.0);

    double intpoints[8][3] = {{-half_Sqrt_three_inv, -half_Sqrt_three_inv, -half_Sqrt_three_inv},
			      {-half_Sqrt_three_inv, -half_Sqrt_three_inv, half_Sqrt_three_inv},
			      {-half_Sqrt_three_inv, half_Sqrt_three_inv, -half_Sqrt_three_inv},
			      {-half_Sqrt_three_inv, half_Sqrt_three_inv, half_Sqrt_three_inv},
			      {half_Sqrt_three_inv, -half_Sqrt_three_inv, -half_Sqrt_three_inv},
			      {half_Sqrt_three_inv, -half_Sqrt_three_inv, half_Sqrt_three_inv},
			      {half_Sqrt_three_inv, half_Sqrt_three_inv, -half_Sqrt_three_inv},
			      {half_Sqrt_three_inv, half_Sqrt_three_inv, half_Sqrt_three_inv}};

    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
	for (int k=0; k<nz; k++){
	  for (int ip=0; ip<8; ip++) {
	    x0[l][0] = x[l][0] = solidlo[0] + delta*(i+0.5+intpoints[ip][0]);
	    x0[l][1] = x[l][1] = solidlo[1] + delta*(j+0.5+intpoints[ip][1]);
	    if (domain->dimension == 3) x0[l][2] = x[l][2] = solidlo[2] + delta*(k+0.5+intpoints[ip][2]);
	    else x0[l][2] = x[l][2] = 0;
	    // Check if the particle is inside cylinder:
	    if (domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2])==1) l++;
	  }
	}
      }
    }
  } else if ((int) input->parsev(args[3]) == 3) {
    // Berstein elements:

    np *= 27;
    mass_ /= 27.0;
    vol_ /= 27.0;
    // Allocate the space in the vectors for np particles:
    grow(np);

    double a = 0.7746/2;

    double intpoints[27][3] = {{-a, -a, -a},
			       {-a, -a, 0},
			       {-a, -a, a},
			       {-a, 0, -a},
			       {-a, 0, 0},
			       {-a, 0, a},
			       {-a, a, -a},
			       {-a, a, 0},
			       {-a, a, a},
			       {0, -a, -a},
			       {0, -a, 0},
			       {0, -a, a},
			       {0, 0, -a},
			       {0, 0, 0},
			       {0, 0, a},
			       {0, a, -a},
			       {0, a, 0},
			       {0, a, a},
			       {a, -a, -a},
			       {a, -a, 0},
			       {a, -a, a},
			       {a, 0, -a},
			       {a, 0, 0},
			       {a, 0, a},
			       {a, a, -a},
			       {a, a, 0},
			       {a, a, a}};

    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
	for (int k=0; k<nz; k++){
	  for (int ip=0; ip<27; ip++) {
	    x0[l][0] = x[l][0] = solidlo[0] + delta*(i+0.5+intpoints[ip][0]);
	    x0[l][1] = x[l][1] = solidlo[1] + delta*(j+0.5+intpoints[ip][1]);
	    if (domain->dimension == 3) x0[l][2] = x[l][2] = solidlo[2] + delta*(k+0.5+intpoints[ip][2]);
	    else x0[l][2] = x[l][2] = 0;
	    // Check if the particle is inside cylinder:
	    if (domain->regions[iregion]->inside(x0[l][0], x0[l][1], x0[l][2])==1) l++;
	  }
	}
      }
    }
  } else {
    cout << "Error: solid command 4th argument should be 1 or 2, but " << (int) input->parsev(args[3]) << "received.\n";
    exit(1);
  }

  np = l; // Adjust np to account for the particles outside the domain
  cout << "np="<< np << endl;

  for (int i=0; i<np;i++) {
    a[i].setZero();
    v[i].setZero();
    f[i].setZero();
    b[i].setZero();
    v_update[i].setZero();
    vol0[i] = vol[i] = vol_;
    rho0[i] = rho[i] = mat->rho0;
    mass[i] = mass_;
    eff_plastic_strain[i] = 0;
    eff_plastic_strain_rate[i] = 0;
    damage[i] = 0;
    damage_init[i] = 0;
    sigma[i].setZero();
    PK1[i].setZero();
    L[i].setZero();
    F[i].setIdentity();
    R[i].setIdentity();
    U[i].setZero();
    D[i].setZero();
    Finv[i].setZero();
    Fdot[i].setZero();

    J[i] = 1;
  }

  if (l!=np) {
    cout << "Error l=" << l << " != np=" << np << endl;
    exit(1);
  }
}


SolCylinder::~SolCylinder()
{

}
