/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

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

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  // cout << "Creating new grid" << endl;

  nx = ny = nz = 0;
  nx_global = ny_global = nz_global = 0;

  cellsize = 0;
  nnodes = 0;

  // Create MPI type for struct Point:
  Point dummy;

  MPI_Datatype type[4] = {MPI_INT, MPI_MPM_TAGINT, MPI_DOUBLE, MPI_INT};
  int blocklen[4] = {1, 1, 3, 3};
  MPI_Aint disp[4];
  MPI_Get_address( &dummy.owner, &disp[0] );
  MPI_Get_address( &dummy.tag, &disp[1] );
  MPI_Get_address( &dummy.x, &disp[2] );
  MPI_Get_address( &dummy.ntype, &disp[3] );
  disp[3] = disp[3] - disp[0];
  disp[2] = disp[2] - disp[0];
  disp[1] = disp[1] - disp[0];
  disp[0] = 0;
  MPI_Type_create_struct(4, blocklen, disp, type, &Pointtype);
  MPI_Type_commit(&Pointtype);

}

Grid::~Grid() {
  // Destroy MPI type:
  MPI_Type_free(&Pointtype);
}

void Grid::init(double *solidlo, double *solidhi) {

#ifdef DEBUG
  std::vector<double> x2plot, y2plot;
#endif

  cout << "In Grid::init() with proc " << universe->me << "\n";
  bool cubic = false;
  bool bernstein = false;
  bool quadratic = false;
  double h = cellsize;

  if (update->method_shape_function.compare("cubic-spline") == 0)
    cubic = true;
  if (update->method_shape_function.compare("quadratic-spline") == 0)
    quadratic = true;
  if (update->method_shape_function.compare("Bernstein-quadratic") == 0) {
    bernstein = true;
    h /= 2;
  }

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  double *boundlo, *boundhi;
  if (update->method->is_TL) {
    boundlo = solidlo;
    boundhi = solidhi;
  } else {
    boundlo = domain->boxlo;
    boundhi = domain->boxhi;
  }

  double Loffsetlo[3] = {MAX(0.0, sublo[0] - boundlo[0]),
                         MAX(0.0, sublo[1] - boundlo[1]),
                         MAX(0.0, sublo[2] - boundlo[2])};
  double Loffsethi_[3] = {MAX(0.0, MIN(subhi[0], boundhi[0]) - boundlo[0]),
                          MAX(0.0, MIN(subhi[1], boundhi[1]) - boundlo[1]),
                          MAX(0.0, MIN(subhi[2], boundhi[2]) - boundlo[2])};

  int noffsetlo[3] = {(int)ceil(Loffsetlo[0] / h), (int)ceil(Loffsetlo[1] / h),
                      (int)ceil(Loffsetlo[2] / h)};

  // cout << "1--- proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";

  int noffsethi_[3] = {(int)ceil(Loffsethi_[0] / h),
                       (int)ceil(Loffsethi_[1] / h),
                       (int)ceil(Loffsethi_[2] / h)};

  // cout << "1--- proc " << universe->me << " noffsethi_=[" << noffsethi_[0] << "," << noffsethi_[1] << "," << noffsethi_[2] << "]\n";

  if (universe->procneigh[0][0] >= 0 &&
      abs(boundlo[0] + noffsetlo[0] * h - sublo[0]) < 1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower x boundary
    // they should belong to procneigh[0][0]
    noffsetlo[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][0] >= 0 &&
      abs(boundlo[1] + noffsetlo[1] * h - sublo[1]) < 1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower y boundary
    // they should belong to procneigh[1][0]
    noffsetlo[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][0] >= 0 &&
      abs(boundlo[2] + noffsetlo[2] * h - sublo[2]) < 1.0e-12) {
    // Some nodes would fall exactly on the subdomain lower z boundary
    // they should belong to procneigh[2][0]
    noffsetlo[2]++;
  }

  if (universe->procneigh[0][1] >= 0 &&
      abs(boundlo[0] + noffsethi_[0] * h - subhi[0]) < 1.0e-12) {
    noffsethi_[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][1] >= 0 &&
      abs(boundlo[1] + noffsethi_[1] * h - subhi[1]) < 1.0e-12) {
    noffsethi_[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][1] >= 0 &&
      abs(boundlo[2] + noffsethi_[2] * h - subhi[2]) < 1.0e-12) {
    noffsethi_[2]++;
  }
  // cout << "2--- proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
  // cout << "2--- proc " << universe->me << " noffsethi_=[" << noffsethi_[0] << "," << noffsethi_[1] << "," << noffsethi_[2] << "]\n";

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

  nx = noffsethi_[0] - noffsetlo[0];
  if (domain->dimension >= 2) {
    ny = noffsethi_[1] - noffsetlo[1];
  } else {
    ny = 1;
  }
  if (domain->dimension >= 3) {
    nz = noffsethi_[2] - noffsetlo[2];
  } else {
    nz = 1;
  }

  if (universe->procneigh[0][1] == -1) {
    while (boundlo[0] + h * (noffsetlo[0] + nx - 0.5) < MIN(subhi[0], boundhi[0]))
      nx++;
  }
  if (universe->procneigh[1][1] == -1) {
    while (boundlo[1] + h * (noffsetlo[1] + ny - 0.5) < MIN(subhi[1], boundhi[1]))
      ny++;
  }
  if (universe->procneigh[2][1] == -1) {
    while (boundlo[2] + h * (noffsetlo[2] + nz - 0.5) < MIN(subhi[2], boundhi[2]))
      nz++;
  }

  // Create nodes that are inside the local subdomain:
  nnodes_local = nx * ny * nz;
  if (nnodes_local < 0) {
    error->one(FLERR,
               "Bad domain decomposition, some CPUs do not have any grid "
               "attached to.\nThis is likely to happen when using TLMPM.\n");
  }
  grow(nnodes_local);

  // #ifdef DEBUG
  cout << "proc " << universe->me << " nx=" << nx << "\tny=" << ny << "\tnz=" << nz <<endl;
  cout << "proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
  cout << "proc " << universe->me << " noffsethi_=[" << noffsethi_[0] << "," << noffsethi_[1] << "," << noffsethi_[2] << "]\n";
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
	} else if (cubic || quadratic) {
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
	pH[l] = 0;
	rigid[l] = false;
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
	plt::annotate(to_string(ntag[l]), x0[l][2], x0[l][1]);
	x2plot.push_back(x0[l][2]);
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
  if (cubic || quadratic || bernstein) delta = 2*h - 1.0e-12;
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
      shared.push_back(ntag[in]);
    }
  }

  nshared = ns.size();
#ifdef DEBUG
  cout << "proc " << universe->me << " has " << ns.size() << " node to send.\n";
#endif

  // Send over the tag and coordinates of the nodes to send:
  int size_n;
  vector<Point> tmp_shared;

  for (int sproc = 0; sproc < universe->nprocs; sproc++) {
    if (sproc == universe->me) {
      size_n = ns.size();
    } else {
      size_n = 0;
    }

    MPI_Bcast(&size_n, 1, MPI_INT, sproc, universe->uworld);

    if (sproc == universe->me) {
      tmp_shared = ns;
    } else {
      tmp_shared.resize(size_n);
    }

    MPI_Bcast(tmp_shared.data(), size_n, Pointtype, sproc, universe->uworld);


    if (sproc == universe->me) {

      // Receive destination lists:

      int size_origin;
      for (int rproc = 0; rproc < universe->nprocs; rproc++) {
        if (rproc != universe->me) {
          MPI_Recv(&size_origin, 1, MPI_INT, rproc, 0, universe->uworld,
                   MPI_STATUS_IGNORE);

          if (size_origin != 0) {
            tagint buf_recv[size_origin];

            MPI_Recv(&buf_recv[0], size_origin, MPI_MPM_TAGINT, rproc, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < size_origin; i++) {
              dest_nshared[rproc].push_back(buf_recv[i]);
            }
          }
        }
      }
    } else {
      vector<tagint> origin_gnodes;

      // Check if the nodes received are in the subdomain:
      for (int i_recv = 0; i_recv < size_n; i_recv++) {
        if (domain->inside_subdomain_extended(tmp_shared[i_recv].x[0],
                                              tmp_shared[i_recv].x[1],
                                              tmp_shared[i_recv].x[2], delta)) {
          gnodes.push_back(tmp_shared[i_recv]);
          origin_gnodes.push_back(tmp_shared[i_recv].tag);
#ifdef DEBUG
          cout << "proc " << universe->me << " received node "
               << tmp_shared[i_recv].tag << " from proc " << sproc << ".\n";
#endif
        }
      }

      // Send origin list:
      int size_origin = origin_gnodes.size();
      MPI_Send(&size_origin, 1, MPI_INT, sproc, 0, universe->uworld);

      if (size_origin != 0) {
        MPI_Send(origin_gnodes.data(), size_origin, MPI_MPM_TAGINT, sproc, 0,
                 MPI_COMM_WORLD);
	origin_nshared[sproc] = origin_gnodes;
      }
    }
  }

#ifdef DEBUG
  // Check that the dest_nshared map is correct:
  for (auto it = dest_nshared.cbegin(); it != dest_nshared.cend(); ++it) {
    cout << "On proc " << universe->me << " proc " << it->first
         << " should receive information of nodes ";
    for (auto ii : it->second) {
      cout << ii << ", ";
    }
    cout << endl;
  }
#endif

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
      cout << "x0[" << ntag[i] << "]=[" << gnodes[in].x[0] << "," << gnodes[in].x[1] << "," << gnodes[in].x[2] << "]\t exits in x=[" << x0[map_ntag[i]][0] << "," << x0[map_ntag[i]][1] << "," << x0[map_ntag[i]][2] << "]\n";
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
    pH[i] = 0;

#ifdef DEBUG
    x2plot.push_back(x0[i][2]);
    y2plot.push_back(x0[i][1]);
    plt::annotate(to_string(ntag[i]), x0[i][2], x0[i][1]);
#endif
  }

#ifdef DEBUG
  plt::plot(x2plot, y2plot, "g+");
#endif
  
  // Determine the total number of nodes:
  bigint nnodes_temp = nnodes_local;
  MPI_Allreduce(&nnodes_temp, &nnodes, 1, MPI_MPM_BIGINT, MPI_SUM, universe->uworld);
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  cout << "Set grid cellsize to " << cellsize << endl;
}

void Grid::grow(int nn){
  //nnodes_local = nn;

  ntag.resize(nn);
  nowner.resize(nn);

  x0.resize(nn);
  x.resize(nn);
  v.resize(nn);
  v_update.resize(nn);
  mb.resize(nn);
  f.resize(nn);
  mass.resize(nn);
  mask.resize(nn);
  pH.resize(nn);

  for (int i=0; i<nn; i++) mask[i] = 1;

  ntype.resize(nn);
  rigid.resize(nn);

}

void Grid::update_grid_velocities()
{
  Vector3d vtemp;
  vtemp.setZero();

  // Update all particles (even the ghost to not have to communicate the result)
  for (int i=0; i<nnodes_local + nnodes_ghost; i++){
    if (!rigid[i]) {
      if (mass[i] != 0) v_update[i] = v[i] + update->dt * (f[i] + mb[i])/mass[i];
      else v_update[i] = v[i];
    } else {
      vtemp = v_update[i]; // v[i] at t-dt
      v_update[i] = v[i];
      v[i] = vtemp;
    }
    // if (update->ntimestep>450)
    // if (isnan(v_update[i](0)))
    //   cout << "update_grid_velocities: in=" << i << ", vn=[" << v[i][0] << "," << v[i][1] << "," << v[i][2] << "], f=[" << f[i][0] << "," << f[i][1] << "," << f[i][2] << "], mb=[" << mb[i][0] << "," << mb[i][1] << "," << mb[i][2] << "], dt=" << update->dt << ", mass[i]=" << mass[i] << endl;
  }
}

void Grid::update_grid_positions()
{
  // Update all particles (even the ghost to not have to communicate the result)
  for (int i=0; i<nnodes_local + nnodes_ghost; i++){
    x[i] += update->dt*v[i];
  }
}

void Grid::reduce_mass_ghost_nodes() {
  vector<double> tmp_mass;
  int j, size_r, size_s, jproc;

  // Some ghost nodes' mass on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list of masses from jproc:

        jproc = idest->first;

	// cout << "proc " << universe->me << " receives size from " << jproc << endl;
        // MPI_Recv(&size_r, 1, MPI_INT, jproc, 0, universe->uworld,
        //          MPI_STATUS_IGNORE);

        // if (size_r != idest->second.size()) {
        //   error->one(
        //       FLERR,
        //       "proc " + to_string(iproc) +
        //           " received a conflicting number of shared nodes from proc " +
        //           to_string(jproc) + "\n");
        // }

	size_r = idest->second.size();

        double buf_recv[size_r];

	// cout << "proc " << universe->me << " receives masses from " << jproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_DOUBLE, jproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Add the received masses to that of the nodes:

        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            mass[map_ntag[j]] += buf_recv[is];
            // cout << "mass[" << j << "]=" << mass[map_ntag[j]] << "\n";
          }
        }
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // me should send the list of masses of all the ghost nodes gotten from
        // iproc to iproc.
        size_s = origin_nshared[iproc].size();

	// cout << "proc " << universe->me << " sends size to " << iproc << endl;
        // MPI_Send(&size_s, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD);

        // Create the list of masses:

        tmp_mass.assign(size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            tmp_mass[is] = mass[map_ntag[j]];
          }
        }

	// cout << "proc " << universe->me << " sends mass to " << iproc << endl;
        MPI_Send(tmp_mass.data(), size_s, MPI_DOUBLE, iproc, 0, MPI_COMM_WORLD);
      }
    }
  }

  // Share the updated masses:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of masses to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

	// cout << "proc " << universe->me << " sends size to " << jproc << endl;
        // MPI_Send(&size_s, 1, MPI_INT, jproc, 0, universe->uworld);

        // Create the list of masses:

        tmp_mass.assign(size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            tmp_mass[is] = mass[map_ntag[j]];
          }
        }

	// cout << "proc " << universe->me << " sends masses to " << jproc << endl;
        MPI_Send(tmp_mass.data(), size_s, MPI_DOUBLE, jproc, 0, MPI_COMM_WORLD);
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

	// cout << "proc " << universe->me << " receives size from " << iproc << endl;
        // MPI_Recv(&size_r, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);

	size_r = origin_nshared[iproc].size();
        double buf_recv[size_r];

	// cout << "proc " << universe->me << " receives masses from " << iproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_DOUBLE, iproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Add the received masses to that of the nodes:

        for (int is = 0; is < size_r; is++) {
	  j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            mass[map_ntag[j]] = buf_recv[is];
            // cout << "mass[" << j << "]=" << mass[map_ntag[j]] << "\n";
          }
        }
      }
    }
  }
}

void Grid::reduce_mass_ghost_nodes_old() {
  tagint in, j, ng;
  vector<tagint> tmp_shared;
  vector<double> tmp_mass, tmp_mass_reduced;

  // MPI_Barrier(universe->uworld);

  for (int proc = 0; proc < universe->nprocs; proc++) {
    if (proc == universe->me) {
      ng = nshared;
    } else {
      ng = 0;
    }

    MPI_Bcast(&ng, 1, MPI_INT, proc, universe->uworld);

    if (proc == universe->me) {
      tmp_shared = shared;
    } else {
      tmp_shared.resize(ng);
    }

    tmp_mass.assign(ng, 0);
    tmp_mass_reduced.assign(ng, 0);

    MPI_Bcast(tmp_shared.data(), ng, MPI_MPM_TAGINT, proc, universe->uworld);

    for (int is = 0; is < ng; is++) {
      j = tmp_shared[is];

      if (map_ntag.count(j)) {
        tmp_mass[is] = mass[map_ntag[j]];
      }
    }

    MPI_Allreduce(tmp_mass.data(), tmp_mass_reduced.data(), ng, MPI_DOUBLE,
                  MPI_SUM, universe->uworld);

    for (int is = 0; is < ng; is++) {
      j = tmp_shared[is];

      if (map_ntag.count(j)) {
        mass[map_ntag[j]] = tmp_mass_reduced[is];
        // cout << "mass[" << map_ntag[j] << "]=" << tmp_mass_reduced[is] <<
        // "\n";
      }
    }

    // // cout << "------------------\n";

    // for (int is=0; is<ng; is++) {
    //   if (proc == universe->me) {
    // 	// Position in array of the shared node: shared[is]
    // 	j = shared[is];
    //   }

    //   // MPI_Barrier(universe->uworld);
    //   MPI_Bcast(&j,1,MPI_MPM_TAGINT,proc,universe->uworld);

    //   if (map_ntag.count(j)) in = map_ntag[j];
    //   else in = -1;

    //   if (in >= 0) {
    // 	mass_local = mass[in];
    //   } else {
    // 	mass_local = 0;
    //   }

    //   MPI_Allreduce(&mass_local, &mass_reduced, 1, MPI_DOUBLE, MPI_SUM,
    //   universe->uworld); if (in >= 0) {
    // 	//mass[in] = mass_reduced;
    // 	cout << "mass[" << in << "]=" << mass_reduced << endl;
    //   }
    // }
  }
}

void Grid::reduce_ghost_nodes(bool only_v) {
  vector<double> tmp;
  int j, k, m, size_r, size_s, jproc, nsend;

  if (only_v)
    nsend = 1 * 3;
  else
    nsend = 3 * 3;

  // Some ghost nodes' mass on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list from jproc:

        jproc = idest->first;

        size_r = idest->second.size();

        double buf_recv[nsend * size_r];

        // cout << "proc " << universe->me << " receives masses from " << jproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_DOUBLE, jproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Add the received data to that of the nodes:

        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            v[m][0] += buf_recv[k];
            v[m][1] += buf_recv[k + 1];
            v[m][2] += buf_recv[k + 2];
            if (!only_v) {
              f[m][0] += buf_recv[k + 3];
              f[m][1] += buf_recv[k + 4];
              f[m][2] += buf_recv[k + 5];

              mb[m][0] += buf_recv[k + 6];
              mb[m][1] += buf_recv[k + 7];
              mb[m][2] += buf_recv[k + 8];
            }
          }
        }
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // me should send the list of data of all the ghost nodes gotten from
        // iproc to iproc.
        size_s = origin_nshared[iproc].size();

        // cout << "proc " << universe->me << " sends size to " << iproc <<
        // endl; MPI_Send(&size_s, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD);

        // Create the list of data:

        tmp.assign(nsend * size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            tmp[k] = v[m][0];
            tmp[k + 1] = v[m][1];
            tmp[k + 2] = v[m][2];
            if (!only_v) {
              tmp[k + 3] = f[m][0];
              tmp[k + 4] = f[m][1];
              tmp[k + 5] = f[m][2];

              tmp[k + 6] = mb[m][0];
              tmp[k + 7] = mb[m][1];
              tmp[k + 8] = mb[m][2];
            }
          }
        }

        // cout << "proc " << universe->me << " sends mass to " << iproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_DOUBLE, iproc, 0,
                 MPI_COMM_WORLD);
      }
    }
  }

  // Share the updated data:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of data to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

        // cout << "proc " << universe->me << " sends size to " << jproc <<
        // endl; MPI_Send(&size_s, 1, MPI_INT, jproc, 0, universe->uworld);

        // Create the list of data:

        tmp.assign(nsend * size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            tmp[k] = v[m][0];
            tmp[k + 1] = v[m][1];
            tmp[k + 2] = v[m][2];
            if (!only_v) {
              tmp[k + 3] = f[m][0];
              tmp[k + 4] = f[m][1];
              tmp[k + 5] = f[m][2];

              tmp[k + 6] = mb[m][0];
              tmp[k + 7] = mb[m][1];
              tmp[k + 8] = mb[m][2];
            }
          }
        }

        // cout << "proc " << universe->me << " sends masses to " << jproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_DOUBLE, jproc, 0,
                 MPI_COMM_WORLD);
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

        // cout << "proc " << universe->me << " receives size from " << iproc <<
        // endl; MPI_Recv(&size_r, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);

        size_r = origin_nshared[iproc].size();
        double buf_recv[nsend * size_r];

        // cout << "proc " << universe->me << " receives masses from " << iproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_DOUBLE, iproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Update the local data with the received values

        for (int is = 0; is < size_r; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            v[m][0] = buf_recv[k];
            v[m][1] = buf_recv[k + 1];
            v[m][2] = buf_recv[k + 2];
            if (!only_v) {
              f[m][0] = buf_recv[k + 3];
              f[m][1] = buf_recv[k + 4];
              f[m][2] = buf_recv[k + 5];

              mb[m][0] = buf_recv[k + 6];
              mb[m][1] = buf_recv[k + 7];
              mb[m][2] = buf_recv[k + 8];
            }
          }
        }
      }
    }
  }
}

void Grid::reduce_ghost_nodes_old(bool only_v) {
  tagint in, j, ng;
  vector<tagint> tmp_shared;
  vector<double> buf, buf_reduced;
  int nsend, k, m;

  if (only_v)
    nsend = 1 * 3;
  else
    nsend = 3 * 3;

  // MPI_Barrier(universe->uworld);

  for (int proc = 0; proc < universe->nprocs; proc++) {
    if (proc == universe->me) {
      ng = nshared;
    } else {
      ng = 0;
    }

    MPI_Bcast(&ng, 1, MPI_INT, proc, universe->uworld);

    if (proc == universe->me) {
      tmp_shared = shared;
    } else {
      tmp_shared.resize(ng);
    }

    MPI_Bcast(tmp_shared.data(), ng, MPI_MPM_TAGINT, proc, universe->uworld);

    buf.assign(nsend * ng, 0);
    buf_reduced.assign(nsend * ng, 0);

    for (int is = 0; is < ng; is++) {
      j = tmp_shared[is];

      if (map_ntag.count(j)) {
        m = map_ntag[j];
        k = nsend * is;
        buf[k] = v[m][0];
        buf[k + 1] = v[m][1];
        buf[k + 2] = v[m][2];
        if (!only_v) {
          buf[k + 3] = f[m][0];
          buf[k + 4] = f[m][1];
          buf[k + 5] = f[m][2];

          buf[k + 6] = mb[m][0];
          buf[k + 7] = mb[m][1];
          buf[k + 8] = mb[m][2];
        }
      }
    }

    MPI_Allreduce(buf.data(), buf_reduced.data(), nsend * ng, MPI_DOUBLE,
                  MPI_SUM, universe->uworld);

    for (int is = 0; is < ng; is++) {
      j = tmp_shared[is];

      if (map_ntag.count(j)) {
        m = map_ntag[j];
        k = nsend * is;
        v[m][0] = buf_reduced[k];
        v[m][1] = buf_reduced[k + 1];
        v[m][2] = buf_reduced[k + 2];
        if (!only_v) {
          f[m][0] = buf_reduced[k + 3];
          f[m][1] = buf_reduced[k + 4];
          f[m][2] = buf_reduced[k + 5];

          mb[m][0] = buf_reduced[k + 6];
          mb[m][1] = buf_reduced[k + 7];
          mb[m][2] = buf_reduced[k + 8];
        }
      }
    }
  }
}

void Grid::reduce_regularized_variables() {
  vector<double> tmp;
  int j, k, m, size_r, size_s, jproc, nsend;

  nsend = 1; // pH

  // Some ghost nodes' mass on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list from jproc:

        jproc = idest->first;

        size_r = idest->second.size();

        double buf_recv[nsend * size_r];

        // cout << "proc " << universe->me << " receives masses from " << jproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_DOUBLE, jproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Add the received data to that of the nodes:

        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            pH[m] += buf_recv[k];
          }
        }
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // me should send the list of data of all the ghost nodes gotten from
        // iproc to iproc.
        size_s = origin_nshared[iproc].size();

        // cout << "proc " << universe->me << " sends size to " << iproc <<
        // endl; MPI_Send(&size_s, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD);

        // Create the list of data:

        tmp.assign(nsend * size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            tmp[k] = pH[m];
          }
        }

        // cout << "proc " << universe->me << " sends mass to " << iproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_DOUBLE, iproc, 0,
                 MPI_COMM_WORLD);
      }
    }
  }

  // Share the updated data:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of data to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

        // cout << "proc " << universe->me << " sends size to " << jproc <<
        // endl; MPI_Send(&size_s, 1, MPI_INT, jproc, 0, universe->uworld);

        // Create the list of data:

        tmp.assign(nsend * size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            tmp[k] = pH[m];
          }
        }

        // cout << "proc " << universe->me << " sends masses to " << jproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_DOUBLE, jproc, 0,
                 MPI_COMM_WORLD);
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

        // cout << "proc " << universe->me << " receives size from " << iproc <<
        // endl; MPI_Recv(&size_r, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);

        size_r = origin_nshared[iproc].size();
        double buf_recv[nsend * size_r];

        // cout << "proc " << universe->me << " receives masses from " << iproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_DOUBLE, iproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Update the local data with the received values

        for (int is = 0; is < size_r; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag.count(j)) {
            m = map_ntag[j];
            k = nsend * is;
            pH[m] = buf_recv[k];
          }
        }
      }
    }
  }
}
