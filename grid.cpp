#include <vector>
#include <math.h>
#include <Eigen/Eigen>
#include <string>
#include "mpm.h"
#include "grid.h"
#include "material.h"
#include "input.h"
#include "memory.h"
#include "update.h"
#include "var.h"
#include "domain.h"
#include "method.h"
#include "universe.h"
#include "error.h"


#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

using namespace std;
using namespace Eigen;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  cout << "Creating new grid" << endl;

  nx = ny = nz = 0;
  nx_global = ny_global = nz_global = 0;

  cellsize = 0;
  nnodes = 0;

  // Create MPI type for struct Point:
  Point dummy;

  MPI_Datatype type[4] = {MPI_INT, MPI_MPM_TAGINT, MPI_DOUBLE, MPI_INT};
  int blocklen[4] = {1, 1, 3, 3};
  MPI_Aint disp[4];
  MPI_Address( &dummy.owner, &disp[0] );
  MPI_Address( &dummy.tag, &disp[1] );
  MPI_Address( &dummy.x, &disp[2] );
  MPI_Address( &dummy.ntype, &disp[3] );
  disp[3] = disp[3] - disp[0];
  disp[2] = disp[2] - disp[0];
  disp[1] = disp[1] - disp[0];
  disp[0] = 0;
  MPI_Type_struct(4, blocklen, disp, type, &Pointtype);
  MPI_Type_commit(&Pointtype);

}

Grid::~Grid()
{
  // Destroy MPI type:
  MPI_Type_free(&Pointtype);
}

void Grid::init(double *solidlo, double *solidhi){

  cout << "In Grid::init()\n";
  bool cubic = false;
  bool bernstein = false;
  double h = cellsize;

  if (update->method_shape_function.compare("cubic-spline")==0) cubic = true;
  if (update->method_shape_function.compare("Bernstein-quadratic")==0) {
    bernstein = true;
    h /= 2;
  }

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


 double Loffsetlo[3] = {MAX(0.0, sublo[0] - boundlo[0]),
			 MAX(0.0, sublo[1] - boundlo[1]),
			 MAX(0.0, sublo[2] - boundlo[2])};

  double Loffsethi[3] = {MIN(0.0, subhi[0] - boundhi[0]),
			 MIN(0.0, subhi[1] - boundhi[1]),
			 MIN(0.0, subhi[2] - boundhi[2])};

  if (Loffsethi[0] > -1.0e-12) Loffsethi[0] = 0;
  if (Loffsethi[1] > -1.0e-12) Loffsethi[1] = 0;
  if (Loffsethi[2] > -1.0e-12) Loffsethi[2] = 0;

  int noffsetlo[3] = {(int) ceil(Loffsetlo[0]/h),
		      (int) ceil(Loffsetlo[1]/h),
		      (int) ceil(Loffsetlo[2]/h)};

  int noffsethi[3] = {(int) ceil(-Loffsethi[0]/h),
		      (int) ceil(-Loffsethi[1]/h),
		      (int) ceil(-Loffsethi[2]/h)};

  if (universe->procneigh[0][0] >= 0 && abs(boundlo[0]+ noffsetlo[0]*h - sublo[0])<1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower x boundary
    // they should below to procneigh[0][0]
    noffsetlo[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][0] >= 0 && abs(boundlo[1]+ noffsetlo[1]*h - sublo[1])<1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower y boundary
    // they should below to procneigh[1][0]
    noffsetlo[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][0] >= 0 && abs(boundlo[2]+ noffsetlo[2]*h - sublo[2])<1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower x boundary
    // they should below to procneigh[2][0]
    noffsetlo[2]++;
  }

  double Lx = (boundhi[0] - boundlo[0]) - (noffsetlo[0] + noffsethi[0])*h;


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

  // Create nodes that are inside the local subdomain:
  nnodes_local = nx*ny*nz;
  grow(nnodes_local);

#ifdef DEBUG
  std::vector<double> x2plot, y2plot;
#endif

  double Lx_global = solidhi[0]-solidlo[0];//+2*cellsize;

  nx_global = ((int) (Lx_global/h))+1;
  while (nx_global*h <= Lx_global+0.5*h) nx_global++;

  if (domain->dimension >= 2) {
    double Ly_global = solidhi[1]-solidlo[1];
    ny_global = ((int) Ly_global/h)+1;
    while (ny_global*h <= Ly_global+0.5*h) ny_global++;
  } else {
    ny_global = 1;
  }

  if (domain->dimension == 3) {
    double Lz_global = solidhi[2]-solidlo[2];
    nz_global = ((int) Lz_global/h)+1;
    while (nz_global*h <= Lz_global+0.5*h) nz_global++;
  } else {
    nz_global = 1;
  }

//   // Need to receive from lower neighbouring CPUs what their nx, ny, and nz:

//   int nx0, ny0, nz0, nx_temp, ny_temp, nz_temp;
//   nx0 = ny0 = nz0 = 0;

//   for (int iproc=0; iproc<universe->nprocs; iproc++){
//     if (iproc == universe->me) {
//       if (universe->procneigh[0][1] >= 0) {
// 	// Send nx to the neighbouring CPU in ascending x
// 	nx_temp = nx0 + nx;
// 	MPI_Send(&nx_temp, 1, MPI_INT, universe->procneigh[0][1], 0, universe->uworld);
//       }

//       if ((domain->dimension >= 2) && (universe->procneigh[1][1] >= 0)){
// 	// Send ny to the neighbouring CPU in ascending y
// 	ny_temp = ny0 + ny;
// 	MPI_Send(&ny_temp, 1, MPI_INT, universe->procneigh[1][1], 0, universe->uworld);
//       }

//       if ((domain->dimension == 3) && (universe->procneigh[2][1] >= 0)){
// 	// Send nz to the neighbouring CPU in ascending z
// 	nz_temp = nz0 + nz;
// 	MPI_Send(&nz_temp, 1, MPI_INT, universe->procneigh[2][1], 0, universe->uworld);
//       }
//     }

//     if (iproc == universe->procneigh[0][0]) {
//       // Receive nx0 from the neighbouring CPU in descending x
//       MPI_Recv(&nx0, 1, MPI_INT, iproc, 0, universe->uworld, MPI_STATUS_IGNORE);
//     }
//     if ((domain->dimension >= 2) && iproc == universe->procneigh[1][0]) {
//       // Receive ny0 from the neighbouring CPU in descending y
//       MPI_Recv(&ny0, 1, MPI_INT, iproc, 0, universe->uworld, MPI_STATUS_IGNORE);
//     }
//     if ((domain->dimension == 3) && iproc == universe->procneigh[2][0]) {
//       // Receive nz0 from the neighbouring CPU in descending y
//       MPI_Recv(&nz0, 1, MPI_INT, iproc, 0, universe->uworld, MPI_STATUS_IGNORE);
//     }
//   }

// #ifdef DEBUG
//   cout << "proc " << universe->me << " nx0=" << nx0 << "\tny0=" << ny0 << "\tnz0=" << nz0 <<endl;
//   cout << "proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
//   if (noffsetlo[0]!=nx0) error->all(FLERR, "noffsetlo[0]!=nx0:" + to_string(noffsetlo[0]) + "!=" + to_string(nx0) + "\n");
//   if (noffsetlo[1]!=ny0) error->all(FLERR, "noffsetlo[1]!=ny0:" + to_string(noffsetlo[1]) + "!=" + to_string(ny0) + "\n");
//   if (noffsetlo[2]!=nz0) error->all(FLERR, "noffsetlo[2]!=nz0:" + to_string(noffsetlo[2]) + "!=" + to_string(nz0) + "\n");
// #endif

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
	} else if (bernstein) {
	  ntype[l][0] = (i+noffsetlo[0]) % 2;
	  ntype[l][1] = (j+noffsetlo[1]) % 2;
	  ntype[l][2] = (k+noffsetlo[2]) % 2;
	} else if (cubic) {
	  ntype[l][0] = min(2,i+noffsetlo[0])-min(nx_global-1-i-noffsetlo[0],2);
	  ntype[l][1] = min(2,j+noffsetlo[1])-min(ny_global-1-j-noffsetlo[1],2);
	  ntype[l][2] = min(2,k+noffsetlo[2])-min(nz_global-1-k-noffsetlo[2],2);
	}

	x[l] = x0[l];
	v[l].setZero();
	v_update[l].setZero();
	f[l].setZero();
	mb[l].setZero();
	mass[l] = 0;
	// R[l].setIdentity();

	ntag[l] = nz_global*ny_global*(i+noffsetlo[0]) + nz_global*(j+noffsetlo[1]) + k+noffsetlo[2];
	// cout << "ntag = " << ntag[l] << endl;

	// Check if ntag[l] already exists:
	if(map_ntag.count(ntag[l]) > 0 ) {
	  error->all(FLERR, "node " + to_string(ntag[l]) + " already exists.");
	}
	map_ntag[ntag[l]] = l;
	nowner[l] = universe->me;

#ifdef DEBUG
	plt::annotate(to_string(ntag[l]), x0[l][0], x0[l][1]);
	x2plot.push_back(x0[l][0]);
	y2plot.push_back(x0[l][1]);
#endif
	l++;
      }
    }
  }

  nnodes_local = l;

  // Give to neighbouring procs ghost nodes:
  vector<Point> ns;
  vector<Point> gnodes;

  double delta;
  if (cubic || bernstein) delta = 2*h - 1.0e-12;
  else delta = h - 1.0e-12;

  cout << "delta=" << delta << endl;

  bool isnt_sublo_boundlo[3] = {true, true, true};
  bool isnt_subhi_boundhi[3] = {true, true, true};

  if (abs(sublo[0] - boundlo[0]) < 1.0e-12) isnt_sublo_boundlo[0] = false;
  if (abs(sublo[1] - boundlo[1]) < 1.0e-12) isnt_sublo_boundlo[1] = false;
  if (abs(sublo[2] - boundlo[2]) < 1.0e-12) isnt_sublo_boundlo[2] = false;

  if (abs(subhi[0] - boundhi[0]) < 1.0e-12) isnt_subhi_boundhi[0] = false;
  if (abs(subhi[1] - boundhi[1]) < 1.0e-12) isnt_subhi_boundhi[1] = false;
  if (abs(subhi[2] - boundhi[2]) < 1.0e-12) isnt_subhi_boundhi[2] = false;

  // For each node, check if it needs to be sent to another proc:
  for (int in=0; in<nnodes_local; in++){
    if (isnt_sublo_boundlo[0] && (x0[in][0] - sublo[0] < delta) ||
	(domain->dimension >= 2) && isnt_sublo_boundlo[1] && (x0[in][1] - sublo[1] < delta) ||
	(domain->dimension == 3) && isnt_sublo_boundlo[2] && (x0[in][2] - sublo[2] < delta) ||
	isnt_subhi_boundhi[0] && (subhi[0] - x0[in][0] < delta) ||
	(domain->dimension >= 2) && isnt_subhi_boundhi[1] && (subhi[1] - x0[in][1] < delta) ||
	(domain->dimension == 3) && isnt_subhi_boundhi[2] && (subhi[2] - x0[in][2] < delta)) {
      Point p = {universe->me, ntag[in], {x0[in][0], x0[in][1], x0[in][2]}, {ntype[in][0], ntype[in][1], ntype[in][2]}};
      ns.push_back(p);
      shared.push_back(in);
    }
  }

  nshared = ns.size();
#ifdef DEBUG
  cout << "proc " << universe->me << " has " << ns.size() << " node to send.\n";
#endif

  MPI_Barrier(universe->uworld);

  // Send over the tag and coordinates of the nodes to send:
  for (int sproc=0; sproc<universe->nprocs; sproc++){
    int size_nr = 0;

    if (sproc == universe->me) {
      int size_ns = ns.size();

      for (int rproc=0; rproc<universe->nprocs; rproc++){
	if (rproc != universe->me) {
	  MPI_Send(&size_ns, 1, MPI_INT, rproc, 0, universe->uworld);
	  MPI_Send(ns.data(), size_ns, Pointtype, rproc, 0, MPI_COMM_WORLD);
	}
      }
    } else {
      MPI_Recv(&size_nr, 1, MPI_INT, sproc, 0, universe->uworld, MPI_STATUS_IGNORE);

      Point buf_recv[size_nr];

      MPI_Recv(&buf_recv[0], size_nr, Pointtype, sproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Check if the nodes received are in the subdomain:
      for (int i_recv=0; i_recv<size_nr; i_recv++) {
	if (domain->inside_subdomain_extended(buf_recv[i_recv].x[0], buf_recv[i_recv].x[1], buf_recv[i_recv].x[2], delta)) {
	  gnodes.push_back(buf_recv[i_recv]);
	}
      }
    }
    MPI_Barrier(universe->uworld);
  }

  nnodes_ghost = gnodes.size();

#ifdef DEBUG
  cout << "proc " << universe->me << " has " << nnodes_ghost << " ghost nodes.\n";
#endif
  grow(nnodes_local + nnodes_ghost);

  // Copy ghost nodes at the end of the local nodes:
  for (int in=0; in<nnodes_ghost; in++) {
    int i = nnodes_local+in;

    nowner[i] = gnodes[in].owner;
    ntag[i] = gnodes[in].tag;

    // Check if ntag[l] already exists:
    if(map_ntag.count(ntag[i]) > 0 ) {
      error->all(FLERR, "node " + to_string(ntag[i]) + " already exists.");
    }

    map_ntag[ntag[i]] = i;

    x0[i][0] = x[i][0] = gnodes[in].x[0];
    x0[i][1] = x[i][1] = gnodes[in].x[1];
    x0[i][2] = x[i][2] = gnodes[in].x[2];

    ntype[i][0] = gnodes[in].ntype[0];
    ntype[i][1] = gnodes[in].ntype[1];
    ntype[i][2] = gnodes[in].ntype[2];

    v[i].setZero();
    v_update[i].setZero();
    f[i].setZero();
    mb[i].setZero();
    mass[i] = 0;

#ifdef DEBUG
    x2plot.push_back(x0[i][0]);
    y2plot.push_back(x0[i][1]);
    plt::annotate(to_string(ntag[i]), x0[i][0], x0[i][1]);
#endif
  }

  MPI_Barrier(universe->uworld);

#ifdef DEBUG
  plt::plot(x2plot, y2plot, "g+");
#endif
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  //nnodes_local = nn;

  string str = "grid-ntag";
  cout << "Growing " << str << endl;
  ntag.resize(nn);

  str = "grid-nowner";
  cout << "Growing " << str << endl;
  nowner.resize(nn);

  str = "grid-x0";
  cout << "Growing " << str << endl;
  x0.resize(nn);

  str = "grid-x";
  cout << "Growing " << str << endl;
  x.resize(nn);

  str = "grid-v";
  cout << "Growing " << str << endl;
  v.resize(nn);

  str = "grid-v_update";
  cout << "Growing " << str << endl;
  v_update.resize(nn);

  str = "grid-mb";
  cout << "Growing " << str << endl;
  mb.resize(nn);

  str = "grid-f";
  cout << "Growing " << str << endl;
  f.resize(nn);

  str = "grid-mass";
  cout << "Growing " << str << endl;
  mass.resize(nn);

  str = "grid-mask";
  cout << "Growing " << str << endl;
  mask.resize(nn);

  for (int i=0; i<nn; i++) mask[i] = 1;

  str = "grid-ntype";
  cout << "Growing " << str << endl;
  ntype.resize(nn);
}

void Grid::update_grid_velocities()
{
  // Update all particles (even the ghost to not have to communicate the result)
  for (int i=0; i<nnodes_local + nnodes_ghost; i++){
    if (mass[i] > 1e-12) v_update[i] = v[i] + update->dt * (f[i] + mb[i])/mass[i];
    else v_update[i] = v[i];
  //   if (ntag[i]==5)
  //   cout << "update_grid_velocities: tag=" << ntag[i] << ", vn=[" << v[i][0] << "," << v[i][1] << "," << v[i][2] << "], f=[" << f[i][0] << "," << f[i][1] << "," << f[i][2] << "], mb=[" << mb[i][0] << "," << mb[i][1] << "," << mb[i][2] << "], dt=" << update->dt << ", mass[i]=" << mass[i] << endl;
  }
}

void Grid::update_grid_positions()
{
  // Update all particles (even the ghost to not have to communicate the result)
  for (int i=0; i<nnodes_local + nnodes_ghost; i++){
    x[i] += update->dt*v[i];
  }
}

void Grid::reduce_mass_ghost_nodes()
{
  tagint in, j, ng;
  double mass_local, mass_reduced;

  MPI_Barrier(universe->uworld);

  for (int proc=0; proc<universe->nprocs; proc++){
    if (proc == universe->me) {
      ng = nshared;
    } else {
      ng = 0;
    }

    MPI_Bcast(&ng,1,MPI_INT,proc,universe->uworld);

    for (int is=0; is<ng; is++) {
      if (proc == universe->me) {
	// Position in array of the shared node: shared[is]
	j = ntag[shared[is]];
      }

      MPI_Barrier(universe->uworld);
      MPI_Bcast(&j,1,MPI_MPM_TAGINT,proc,universe->uworld);

      if (map_ntag.count(j)) in = map_ntag[j];
      else in = -1;

      if (in >= 0) {
	mass_local = mass[in];	
      } else {
	mass_local = 0;
      }

      MPI_Allreduce(&mass_local, &mass_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
      if (in >= 0) mass[in] = mass_reduced;
    }
  }
}

void Grid::reduce_ghost_nodes(bool only_v)
{
  tagint in, j, ng;
  double f_local[3], f_reduced[3], v_local[3], v_reduced[3];

  MPI_Barrier(universe->uworld);

  for (int proc=0; proc<universe->nprocs; proc++){
    if (proc == universe->me) {
      ng = nshared;
    } else {
      ng = 0;
    }

    MPI_Bcast(&ng,1,MPI_INT,proc,universe->uworld);

    for (int is=0; is<ng; is++) {
      if (proc == universe->me) {
	// Position in array of the shared node: shared[is]
	j = ntag[shared[is]];
      }

      MPI_Barrier(universe->uworld);
      MPI_Bcast(&j,1,MPI_MPM_TAGINT,proc,universe->uworld);

      if (map_ntag.count(j)) in = map_ntag[j];
      else in = -1;

      if (in >= 0) {
	if (!only_v) {
	  f_local[0] = f[in][0];
	  f_local[1] = f[in][1];
	  f_local[2] = f[in][2];
	}
	v_local[0] = v[in][0];
	v_local[1] = v[in][1];
	v_local[2] = v[in][2];	
      } else {
	if (!only_v) {
	  f_local[0] = f_local[1] = f_local[2] = 0;
	}
	v_local[0] = v_local[1] = v_local[2] = 0;
      }

      if (!only_v) {
	MPI_Allreduce(f_local, f_reduced, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);

	if (in >= 0) {
	  f[in][0] = f_reduced[0];
	  f[in][1] = f_reduced[1];
	  f[in][2] = f_reduced[2];
	}
      }

      MPI_Allreduce(v_local, v_reduced, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);

      if (in >= 0) {
	v[in][0] = v_reduced[0];
	v[in][1] = v_reduced[1];
	v[in][2] = v_reduced[2];
      }
    }
  }
}
