#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include <vector>
#include "memory.h"
#include "update.h"
#include "var.h"
#include "domain.h"
#include "method.h"
#include "error.h"

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;

  x = x0 = NULL;
  v = v_update = NULL;
  mb = f = NULL;

  mass = NULL;
  mask = NULL;
  ntype = NULL;

  // R = NULL;

  // C = NULL;

  cellsize = 0;
  nnodes = 0;
}

Grid::~Grid()
{
  memory->destroy(x);
  memory->destroy(x0);
  memory->destroy(v);
  memory->destroy(v_update);
  memory->destroy(mb);
  memory->destroy(f);
  memory->destroy(mass);
  memory->destroy(mask);
  memory->destroy(ntype);
}

void Grid::init(double *solidlo, double *solidhi){

  double h = cellsize;

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  double *boundlo, *boundhi;
  if (update->method->is_TL) {
    boundlo = solidlo;
    boundhi = solidhi;
  }
  else {
    boundlo = domain->boxlo;
    boundhi = domain->boxhi;
  }

  bool is_cubic = false;

  if (update->method_shape_function.compare("cubic-spline")==0) is_cubic = true;

  double Loffsetlo[3] = {MAX(0.0, domain->sublo[0] - boundlo[0] - is_cubic*2*h),
			 MAX(0.0, domain->sublo[1] - boundlo[1] - is_cubic*2*h),
			 MAX(0.0, domain->sublo[2] - boundlo[2] - is_cubic*2*h)};

  double Loffsethi[3] = {MIN(0.0, domain->subhi[0] - boundhi[0] + is_cubic*2*h),
			 MIN(0.0, domain->subhi[1] - boundhi[1] + is_cubic*2*h),
			 MIN(0.0, domain->subhi[2] - boundhi[2] + is_cubic*2*h)};

  int noffsetlo[3] = {(int) (Loffsetlo[0]/h),
		      (int) (Loffsetlo[1]/h),
		      (int) (Loffsetlo[2]/h)};

  int noffsethi[3] = {(int) (-Loffsethi[0]/h),
		      (int) (-Loffsethi[1]/h),
		      (int) (-Loffsethi[2]/h)};

  double Lx = (boundhi[0] - boundlo[0]) - (noffsetlo[0] + noffsethi[0])*h;

  if (update->method_shape_function.compare("Bernstein-quadratic")==0)
    h /= 2;

  nx = ((int) (Lx/h))+1;
  while (nx*h <= Lx+0.5*h) nx++;

  if (domain->dimension >= 2) {
    double Ly = (boundhi[1] - boundlo[1]) - (noffsetlo[1] + noffsethi[1])*h;
    ny = ((int) Ly/h)+1;
    while (ny*h <= Ly+0.5*h) ny++;   
   } else {
    ny = 1;
   }

  if (domain->dimension == 3) {
    double Lz = (boundhi[2] - boundlo[2]) - (noffsetlo[2] + noffsethi[2])*h;
    nz = ((int) Lz/h)+1;
    while (nz*h <= Lz+0.5*h) nz++;   
   } else {
    nz = 1;
   }


  int nn = nx*ny*nz;
  grow(nn);

#ifdef DEBUG
  std::vector<double> x2plot, y2plot;
#endif


  int l=0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
	x0[l][0] = boundlo[0] + (noffsetlo[0] + i)*h;//h*(i-1);
	if (domain->dimension >= 2) x0[l][1] = boundlo[1] + (noffsetlo[1] + j)*h;
	else x0[l][1] = 0;
	if (domain->dimension == 3) x0[l][2] = boundlo[2] + (noffsetlo[2] + k)*h;
	else x0[l][2] = 0;

	if (update->method_shape_function.compare("linear")==0) {
	  ntype[l][0] = 0;
	  ntype[l][1] = 0;
	  ntype[l][2] = 0;
	} else if (update->method_shape_function.compare("Bernstein-quadratic")==0) {
	  ntype[l][0] = i % 2;
	  ntype[l][1] = j % 2;
	  ntype[l][2] = k % 2;
	} else if (update->method_shape_function.compare("cubic-spline")==0) {
	  ntype[l][0] = min(2,i)-min(nx-1-i,2);
	  ntype[l][1] = min(2,j)-min(ny-1-j,2);
	  ntype[l][2] = min(2,k)-min(nz-1-k,2);
	}
	
	x[l] = x0[l];
	v[l].setZero();
	v_update[l].setZero();
	f[l].setZero();
	mb[l].setZero();
	mass[l] = 0;
	// R[l].setIdentity();

#ifdef DEBUG
	x2plot.push_back(x0[l][0]);
	y2plot.push_back(x0[l][1]);
#endif
	l++;
      }
    }
  }

  
#ifdef DEBUG
  plt::plot(x2plot, y2plot, "g+");
#endif
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  nnodes = nn;

  if (x0 == NULL) x0 = new Eigen::Vector3d[nn];

  if (x == NULL) x = new Eigen::Vector3d[nn];

  if (v == NULL) v = new Eigen::Vector3d[nn];

  if (v_update == NULL) v_update = new Eigen::Vector3d[nn];

  if (mb == NULL) mb = new Eigen::Vector3d[nn];

  if (f == NULL) f = new Eigen::Vector3d[nn];

  // if (R == NULL) R = new Eigen::Matrix3d[nn];

  // if (C == NULL) C = new Eigen::Vector3d[nn];

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
    if (mass[i] > 1e-12) v_update[i] = v[i] + update->dt * (f[i] + mb[i])/mass[i];
    else v_update[i] = v[i];
    // if (update->ntimestep>450)
    //   if (i==0)
    // 	cout << "update_grid_velocities: in=" << i << ", vn=[" << v[i][0] << "," << v[i][1] << "," << v[i][2] << "], f=[" << f[i][0] << "," << f[i][1] << "," << f[i][2] << "], b=[" << b[i][0] << "," << b[i][1] << "," << b[i][2] << "], dt=" << update->dt << ", mass[i]=" << mass[i] << endl;
  }
}

void Grid::update_grid_positions()
{
  for (int i=0; i<nnodes; i++){
    x[i] += update->dt*v[i];
  }
}
