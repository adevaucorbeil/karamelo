#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include <vector>
#include "memory.h"
#include "update.h"
#include "var.h"

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;

  x = x0 = NULL;
  v = v_update = NULL;
  b = f = NULL;

  mass = NULL;
  mask = NULL;
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
}

void Grid::init(double *solidlo, double *solidhi){
  double Lx = solidhi[0]-solidlo[0];
  double Ly = solidhi[1]-solidlo[1];
  double Lz = solidhi[2]-solidlo[2];
  int nx = (int) Lx/cellsize;
  int ny = (int) Ly/cellsize;
  int nz = (int) Lz/cellsize;
  while ((nx-2)*cellsize <= Lx + 0.5*cellsize) nx++;
  while ((ny-2)*cellsize <= Ly + 0.5*cellsize) ny++;
  while ((nz-2)*cellsize <= Lz + 0.5*cellsize) nz++;

  int nn = nx*ny*nz;
  grow(nn);


  int l=0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	x0[l][0] = solidlo[0] + cellsize*(i-1);
	x0[l][1] = solidlo[1] + cellsize*(j-1);
	x0[l][2] = solidlo[2] + cellsize*(k-1);

	x[l] = x0[l];
	v[l].setZero();
	v_update[l].setZero();
	f[l].setZero();
	b[l].setZero();
	mass[l] = 0;

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

  string str = "grid-mass";
  cout << "Growing " << str << endl;
  mass = memory->grow(mass, nn, str);

  str = "grid-mask";
  cout << "Growing " << str << endl;
  mask = memory->grow(mask, nn, str);

  for (int i=0; i<nn; i++) mask[i] = 1;
}

void Grid::update_grid_velocities()
{
  for (int i=0; i<nnodes; i++){
    if (mass[i] > 0) v_update[i] = v[i] + update->dt/mass[i] * (f[i] + b[i]);
    else v_update[i] = v[i];
  }
}

void Grid::update_grid_positions()
{
  for (int i=0; i<nnodes; i++){
    x[i] += update->dt*v_update[i];
  }
}
