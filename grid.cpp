#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include <vector>
#include "memory.h"
#include "update.h"
#include "var.h"
#include "domain.h"

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  // Check that a method is available:
  if (update->method == NULL) {
    cout << "Error: a method should be defined before creating a solid!" << endl;
    exit(1);
  }

  cout << "Creating new grid" << endl;

  x = x0 = NULL;
  v = v_update = NULL;
  b = f = NULL;

  mass = NULL;
  mask = NULL;
  ntype = NULL;

  R = NULL;

  C = NULL;
}

Grid::~Grid()
{
  memory->destroy(x);
  memory->destroy(x0);
  memory->destroy(v);
  memory->destroy(v_update);
  memory->destroy(b);
  memory->destroy(f);
  memory->destroy(mass);
  memory->destroy(mask);
  memory->destroy(ntype);
  memory->destroy(R);
  memory->destroy(C);
}

void Grid::init(double *solidlo, double *solidhi){
  double Lx = solidhi[0]-solidlo[0];//+2*cellsize;
  double Ly = solidhi[1]-solidlo[1];//+2*cellsize;

  double h = cellsize;

  if (update->method_shape_function.compare("Bernstein-quadratic")==0)
    h /= 2;

  nx = ((int) Lx/h)+1;
  ny = ((int) Ly/h)+1;

  //while (nx*h <= Lx+0.5*h) nx++;
  //while (ny*h <= Ly+0.5*h) ny++;

  if (domain->dimension == 3) {
    double Lz = solidhi[2]-solidlo[2]+2*h;
    nz = ((int) Lz/h)+1;
    //while (nz*h <= Lz+0.5*h) nz++;   
   } else {
    nz = 1;
   }


  int nn = nx*ny*nz;
  grow(nn);


  int l=0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	x0[l][0] = solidlo[0] + i*h;//h*(i-1);
	x0[l][1] = solidlo[1] + j*h;//h*(j-1);
	if (domain->dimension == 3) x0[l][2] = solidlo[2] + k*h;//h*(k-1);
	else x0[l][2] = 0;

	ntype[l][0] = i % 2;
	ntype[l][1] = j % 2;
	ntype[l][2] = k % 2;

	x[l] = x0[l];
	v[l].setZero();
	v_update[l].setZero();
	f[l].setZero();
	b[l].setZero();
	mass[l] = 0;
	R[l].setIdentity();

	l++;
      }
    }
  }

  
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  nnodes = nn;

  if (x0 == NULL) x0 = new Eigen::Vector3d[nn];
  else {
    cout << "Error: x0 already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (x == NULL) x = new Eigen::Vector3d[nn];
  else {
    cout << "Error: x already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (v == NULL) v = new Eigen::Vector3d[nn];
  else {
    cout << "Error: v already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (v_update == NULL) v_update = new Eigen::Vector3d[nn];
  else {
    cout << "Error: v_update already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (b == NULL) b = new Eigen::Vector3d[nn];
  else {
    cout << "Error: b already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (f == NULL) f = new Eigen::Vector3d[nn];
  else {
    cout << "Error: f already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (R == NULL) R = new Eigen::Matrix3d[nn];
  else {
    cout << "Error: R already exists, I don't know how to grow it!\n";
    exit(1);
  }

  if (C == NULL) C = new Eigen::Vector3d[ncells];
  else {
    cout << "Error: C already exists, I don't know how to grow it!\n";
    exit(1);
  }

  string str = "grid-mass";
  cout << "Growing " << str << endl;
  mass = memory->grow(mass, nn, str);

  str = "grid-mask";
  cout << "Growing " << str << endl;
  mask = memory->grow(mask, nn, str);

  for (int i=0; i<nn; i++) mask[i] = 1;

  str = "grid-ntype";
  cout << "Growing " << str << endl;
  ntype = memory->grow(ntype, nn, 3, str);
}

void Grid::update_grid_velocities()
{
  for (int i=0; i<nnodes; i++){
    if (mass[i] > 0) v_update[i] = v[i] + update->dt * (f[i]/mass[i] + b[i]);
    else v_update[i] = v[i];
    // if (update->ntimestep>450)
    //   if (i==0)
    // 	cout << "update_grid_velocities: in=" << i << ", vn=[" << v[i][0] << "," << v[i][1] << "," << v[i][2] << "], f=[" << f[i][0] << "," << f[i][1] << "," << f[i][2] << "], b=[" << b[i][0] << "," << b[i][1] << "," << b[i][2] << "], dt=" << update->dt << ", mass[i]=" << mass[i] << endl;
  }
}

void Grid::update_grid_positions()
{
  for (int i=0; i<nnodes; i++){
    x[i] += update->dt*v_update[i];
  }
}
