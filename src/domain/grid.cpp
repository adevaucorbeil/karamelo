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
#include <matrix.h>
#include <string>
#include <mpm.h>
#include <grid.h>
#include <material.h>
#include <input.h>
#include <memory.h>
#include <update.h>
#include <var.h>
#include <domain.h>
#include <method.h>
#include <universe.h>
#include <error.h>

using namespace std;


Grid::Grid(MPM *mpm) :
  Pointers(mpm)
{
  // cout << "Creating new grid" << endl;

  nx = ny = nz = 0;
  nx_global = ny_global = nz_global = 0;

  cellsize = 0;
  nnodes = 0;
}

void Grid::init(float *solidlo, float *solidhi) {

  bool linear = false;
  bool cubic = false;
  bool bernstein = false;
  bool quadratic = false;
  float h = cellsize;

  if (update->shape_function == Update::ShapeFunctions::LINEAR)
    linear = true;
  if (update->shape_function == Update::ShapeFunctions::CUBIC_SPLINE)
    cubic = true;
  if (update->shape_function == Update::ShapeFunctions::QUADRATIC_SPLINE)
    quadratic = true;
  if (update->shape_function == Update::ShapeFunctions::BERNSTEIN) {
    bernstein = true;
    h /= 2;
  }

  float *sublo = domain->sublo;
  float *subhi = domain->subhi;

  float *boundlo, *boundhi;
  if (update->method->is_TL) {
    boundlo = solidlo;
    boundhi = solidhi;
  } else {
    boundlo = domain->boxlo;
    boundhi = domain->boxhi;
  }

  float Loffsetlo[3] = {MAX(0.0f, sublo[0] - boundlo[0]),
			MAX(0.0f, sublo[1] - boundlo[1]),
			MAX(0.0f, sublo[2] - boundlo[2])};
  float Loffsethi_[3] = {MAX(0.0f, MIN(subhi[0], boundhi[0]) - boundlo[0]),
			 MAX(0.0f, MIN(subhi[1], boundhi[1]) - boundlo[1]),
			 MAX(0.0f, MIN(subhi[2], boundhi[2]) - boundlo[2])};

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
      abs(boundlo[0] + noffsethi_[0] * h - MIN(subhi[0], boundhi[0])) < 1.0e-12) {
    noffsethi_[0]++;
  }
  if (domain->dimension >= 2 && universe->procneigh[1][1] >= 0 &&
      abs(boundlo[1] + noffsethi_[1] * h - MIN(subhi[1], boundhi[1])) < 1.0e-12) {
    noffsethi_[1]++;
  }
  if (domain->dimension == 3 && universe->procneigh[2][1] >= 0 &&
      abs(boundlo[2] + noffsethi_[2] * h - MIN(subhi[2], boundhi[2])) < 1.0e-12) {
    noffsethi_[2]++;
  }
  // cout << "2--- proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
  // cout << "2--- proc " << universe->me << " noffsethi_=[" << noffsethi_[0] << "," << noffsethi_[1] << "," << noffsethi_[2] << "]\n";

  float Lx_global = solidhi[0]-solidlo[0];//+2*cellsize;

  nx_global = ((int) (Lx_global/h))+1;
  while (nx_global*h <= Lx_global+0.5*h) nx_global++;

  if (domain->dimension >= 2) {
    float Ly_global = solidhi[1]-solidlo[1];
    ny_global = ((int) Ly_global/h)+1;
    while (ny_global*h <= Ly_global+0.5*h) ny_global++;
  } else {
    ny_global = 1;
  }

  if (domain->dimension == 3) {
    float Lz_global = solidhi[2]-solidlo[2];
    nz_global = ((int) Lz_global/h)+1;
    while (nz_global*h <= Lz_global+0.5*h) nz_global++;
  } else {
    nz_global = 1;
  }

  // Determine the total number of nodes:
  nnodes = nx_global * ny_global * nz_global;
  
  map_ntag = Kokkos::View<tagint*>("map_ntag", nnodes);
  Kokkos::View<tagint*> map_ntag = this->map_ntag;
  Kokkos::parallel_for("set map_ntag", nnodes,
  KOKKOS_LAMBDA (const int &in)
  {
    map_ntag[in] = -1;
  });

  nx = MAX(0, noffsethi_[0] - noffsetlo[0]);
  if (domain->dimension >= 2) {
    ny = MAX(0, noffsethi_[1] - noffsetlo[1]);
  } else {
    ny = 1;
  }
  if (domain->dimension >= 3) {
    nz = MAX(0, noffsethi_[2] - noffsetlo[2]);
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

#ifdef DEBUG
  cout << "proc " << universe->me << " nx=" << nx << "\tny=" << ny << "\tnz=" << nz <<endl;
  cout << "proc " << universe->me << " noffsetlo=[" << noffsetlo[0] << "," << noffsetlo[1] << "," << noffsetlo[2] << "]\n";
  cout << "proc " << universe->me << " noffsethi_=[" << noffsethi_[0] << "," << noffsethi_[1] << "," << noffsethi_[2] << "]\n";
#endif

  // Create nodes that are inside the local subdomain:
  nnodes_local = nx * ny * nz;

  if (!update->method->is_TL && nnodes_local <= 0) {
    error->one(FLERR,
	       "Bad domain decomposition, some CPUs do not have any grid "
	       "attached to.\n");
  }

  // Give to neighbouring procs ghost nodes:
  vector<Point> ns;
  vector<Point> gnodes;

  float delta;
  if (cubic || quadratic || bernstein) delta = 2*h - 1.0e-12;
  else delta = h - 1.0e-12;

  if (universe->me == 0) {
    cout << "delta=" << delta << endl;
  }

  bool isnt_sublo_boundlo[3] = {true, true, true};
  bool isnt_subhi_boundhi[3] = {true, true, true};

  if (boundlo[0] - sublo[0] > -1.0e-12) isnt_sublo_boundlo[0] = false;
  if (boundlo[1] - sublo[1] > -1.0e-12) isnt_sublo_boundlo[1] = false;
  if (boundlo[2] - sublo[2] > -1.0e-12) isnt_sublo_boundlo[2] = false;

  if (subhi[0] - boundhi[0] > -1.0e-12) isnt_subhi_boundhi[0] = false;
  if (subhi[1] - boundhi[1] > -1.0e-12) isnt_subhi_boundhi[1] = false;
  if (subhi[2] - boundhi[2] > -1.0e-12) isnt_subhi_boundhi[2] = false;

  // For each node, check if it needs to be sent to another proc:
  Kokkos::View<tagint*>::HostMirror ntag_mirror = create_mirror(ntag);
  Kokkos::deep_copy(ntag_mirror, ntag);

  Kokkos::View<Vector3d*>::HostMirror x0_mirror = create_mirror(x0);
  Kokkos::deep_copy(x0_mirror, x0);

  Kokkos::View<Vector3i*>::HostMirror ntype_mirror = create_mirror(ntype);
  Kokkos::deep_copy(ntype_mirror, ntype);

  for (int in=0; in<nnodes_local; in++){
    if (isnt_sublo_boundlo[0] && (x0_mirror[in][0] - sublo[0] < delta) ||
	(domain->dimension >= 2) && isnt_sublo_boundlo[1] && (x0_mirror[in][1] - sublo[1] < delta) ||
	(domain->dimension == 3) && isnt_sublo_boundlo[2] && (x0_mirror[in][2] - sublo[2] < delta) ||
	isnt_subhi_boundhi[0] && (subhi[0] - x0_mirror[in][0] < delta) ||
	(domain->dimension >= 2) && isnt_subhi_boundhi[1] && (subhi[1] - x0_mirror[in][1] < delta) ||
	(domain->dimension == 3) && isnt_subhi_boundhi[2] && (subhi[2] - x0_mirror[in][2] < delta)) {
      Point p = {universe->me, ntag_mirror[in], {x0_mirror[in][0], x0_mirror[in][1], x0_mirror[in][2]},
      {ntype_mirror[in][0], ntype_mirror[in][1], ntype_mirror[in][2]}};
      ns.push_back(p);
      shared.push_back(ntag_mirror[in]);
    }
  }

  nshared = ns.size();

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

    MPI_Bcast(tmp_shared.data(), size_n, universe->Pointtype, sproc, universe->uworld);


    if (sproc == universe->me) {

      // Receive destination lists:

      int size_origin;
      for (int rproc = 0; rproc < universe->nprocs; rproc++) {
        if (rproc != universe->me) {
          MPI_Recv(&size_origin, 1, MPI_INT, rproc, 0, universe->uworld,
                   MPI_STATUS_IGNORE);

          if (size_origin != 0) {
            vector<tagint> buf_recv(size_origin);

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

  nnodes_ghost = gnodes.size();

  grow(nnodes_local + nnodes_ghost);

  int me = universe->me;
  float boundlo_0   = boundlo[0],   boundlo_1   = boundlo[1],   boundlo_2   = boundlo[2];
  int noffsetlo_0 = noffsetlo[0], noffsetlo_1 = noffsetlo[1], noffsetlo_2 = noffsetlo[2];
  int dimension = domain->dimension;
  int nx = this->nx, ny = this->ny, nz = this->nz;
  int nx_global = this->nx_global, ny_global = this->ny_global, nz_global = this->nz_global;

  Kokkos::View<tagint*> ntag = this->ntag;
  Kokkos::View<int*> nowner = this->nowner;

  Kokkos::View<Vector3d*> x = this->x;
  Kokkos::View<Vector3d*> x0 = this->x0;
  Kokkos::View<Vector3d*> v = this->v;
  Kokkos::View<Vector3d*> v_update = this->v_update;
  Kokkos::View<Vector3d*> mb = this->mb;
  Kokkos::View<Vector3d*> f = this->f;

  Kokkos::View<float*> mass = this->mass;
  Kokkos::View<bool*> rigid = this->rigid;
  Kokkos::View<Vector3i*> ntype = this->ntype;

  Kokkos::parallel_for(__PRETTY_FUNCTION__, Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
    { 0, 0, 0 }, { (size_t)nx, (size_t)ny, (size_t)nz }),
  KOKKOS_LAMBDA (int i, int j, int k)
  {
    int l = i + nx*(j + ny*k);

    x0[l][0] = boundlo_0 + (noffsetlo_0 + i)*h;//h*(i-1);
    if (dimension >= 2) x0[l][1] = boundlo_1 + (noffsetlo_1 + j)*h;
    else x0[l][1] = 0;
    if (dimension == 3) x0[l][2] = boundlo_2 + (noffsetlo_2 + k)*h;
    else x0[l][2] = 0;

    if (linear) {
      ntype[l][0] = 0;
      ntype[l][1] = 0;
      ntype[l][2] = 0;
    } else if (bernstein) {
      ntype[l][0] = (i+noffsetlo_0)%2;
      ntype[l][1] = (j+noffsetlo_1)%2;
      ntype[l][2] = (k+noffsetlo_2)%2;
    } else if (cubic || quadratic) {
      ntype[l][0] = MIN(2,i+noffsetlo_0)-MIN(nx_global-1-i-noffsetlo_0,2);
      ntype[l][1] = MIN(2,j+noffsetlo_1)-MIN(ny_global-1-j-noffsetlo_1,2);
      ntype[l][2] = MIN(2,k+noffsetlo_2)-MIN(nz_global-1-k-noffsetlo_2,2);
    }

    x[l] = x0[l];
    v[l] = Vector3d();
    v_update[l] = Vector3d();
    f[l] = Vector3d();
    mb[l] = Vector3d();
    mass[l] = 0;
    rigid[l] = false;
    // R[l].setIdentity();

    ntag[l] = nz_global*ny_global*(i+noffsetlo_0) + nz_global*(j+noffsetlo_1) + k+noffsetlo_2;
    // cout << "ntag = " << ntag[l] << endl;

    // Check if ntag[l] already exists:
    // if(map_ntag[ntag[l]] != -1 ) {
    //   error->all(FLERR, "node " + to_string(ntag[l]) + " already exists.");
    // }
    map_ntag[ntag[l]] = l;
    nowner[l] = me;
  });

  // Copy ghost nodes at the end of the local nodes:
  for (int in=0; in<nnodes_ghost; in++) {
    int i = nnodes_local+in;

    nowner[i] = gnodes[in].owner;
    ntag[i] = gnodes[in].tag;

    // Check if ntag[l] already exists:
    if(map_ntag[ntag[i]] != -1 ) {
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

    v[i] = Vector3d();
    v_update[i] = Vector3d();
    f[i] = Vector3d();
    mb[i] = Vector3d();
    mass[i] = 0;
  }
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  if (universe->me == 0) {
    cout << "Set grid cellsize to " << cellsize << endl;
  }
}

void Grid::grow(int nn){
  //nnodes_local = nn;

  ntag   = Kokkos::View<tagint*>("ntag",   nn);
  nowner = Kokkos::View<int*>   ("nowner", nn);

  x0       = Kokkos::View<Vector3d*>("x0",       nn);
  x        = Kokkos::View<Vector3d*>("x",        nn);
  v        = Kokkos::View<Vector3d*>("v",        nn);
  v_update = Kokkos::View<Vector3d*>("v_update", nn);
  mb       = Kokkos::View<Vector3d*>("mb",       nn);
  f        = Kokkos::View<Vector3d*>("f",        nn);

  mass = Kokkos::View<float*>("mass", nn);
  mask = Kokkos::View<int*>   ("mask", nn);

  Kokkos::View<int*> mask = this->mask;

  Kokkos::parallel_for(__PRETTY_FUNCTION__, nn,
  KOKKOS_LAMBDA(int i)
  {
    mask[i] = 1;
  });

  ntype = Kokkos::View<Vector3i*>("ntype", nn);
  rigid = Kokkos::View<bool*>    ("rigid", nn);

  T        = Kokkos::View<float*>("T",        nn);
  T_update = Kokkos::View<float*>("T_update", nn);
  Qext     = Kokkos::View<float*>("Qext",     nn);
  Qint     = Kokkos::View<float*>("Qint",     nn);
}


void Grid::reduce_mass_ghost_nodes() {
  vector<vector<float>> tmp_mass_vect(universe->nprocs);
  vector<vector<float>> buf_recv_vect(universe->nprocs);
  //vector<float> tmp_mass;
  int j, size_r, size_s, jproc;

  // 1. Pack mass of nodes to be sent back to their owner:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {
      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        jproc = idest->first;
	size_r = idest->second.size();

        buf_recv_vect[jproc].resize(size_r);
      }
    } else {
      // Check if iproc is listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // me should send the list of masses of all the ghost nodes gotten from
        // iproc to iproc.
        size_s = origin_nshared[iproc].size();

        // Create the list of masses:

        tmp_mass_vect[iproc].resize(size_s);
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            tmp_mass_vect[iproc][is] = mass[map_ntag[j]];
          } else {
	    error->one(FLERR, "Grid node j does not exist on this CPU.\n");
	  }
        }
      }      
    }
  }

  // 2. New send and receive
  for (int i = 0; i < universe->sendnrecv.size(); i++) {
    if (universe->sendnrecv[i][0] == 0) {
      // Receive
      jproc = universe->sendnrecv[i][1];

      MPI_Recv(&buf_recv_vect[jproc][0], dest_nshared[jproc].size(), MPI_FLOAT,
               jproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      // Send
      jproc = universe->sendnrecv[i][1];
      MPI_Send(tmp_mass_vect[jproc].data(), origin_nshared[jproc].size(),
               MPI_FLOAT, jproc, 0, MPI_COMM_WORLD);
    }
  }

  // 3. Add the received masses to that of the nodes:
  for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
       ++idest) {
    for (int is = 0; is < idest->second.size(); is++) {
      j = idest->second[is];
      if (map_ntag[j] != -1) {
	mass[map_ntag[j]] += buf_recv_vect[idest->first][is];
      } else {
	error->one(FLERR, "Grid node j does not exist on this CPU.\n");
      }
    }
  }

  // Share the updated masses:
  // 1. Pack
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of masses to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

	// cout << "proc " << universe->me << " sends size to " << jproc << endl;

        // Create the list of masses:
        tmp_mass_vect[jproc].resize(size_s);
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            tmp_mass_vect[jproc][is] = mass[map_ntag[j]];
          } else {
	    error->one(FLERR, "Grid node j does not exist on this CPU.\n");
	  }
        }
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

	// cout << "proc " << universe->me << " receives size from " << iproc << endl;

	size_r = origin_nshared[iproc].size();
        buf_recv_vect[iproc].resize(size_r);
      }
    }
  }


  // 2. New send and receive
  for (int i = 0; i < universe->sendnrecv.size(); i++) {
    if (universe->sendnrecv[i][0] == 0) {
      // Receive
      jproc = universe->sendnrecv[i][1];

      MPI_Send(tmp_mass_vect[jproc].data(), dest_nshared[jproc].size(), MPI_FLOAT,
               jproc, 0, MPI_COMM_WORLD);
    } else {
      // Send
      jproc = universe->sendnrecv[i][1];
      MPI_Recv(&buf_recv_vect[jproc][0], origin_nshared[jproc].size(),
	       MPI_FLOAT, jproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // 3. Overwrite the nodes' mass with the received values
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc != universe->me) {
      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
	size_r = origin_nshared[iproc].size();
        for (int is = 0; is < size_r; is++) {
	  j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            mass[map_ntag[j]] = buf_recv_vect[iproc][is];
          }
        }
      }
    }
  }
}


void Grid::reduce_mass_ghost_nodes_old() {
  vector<float> tmp_mass;
  int j, size_r, size_s, jproc;

  // Some ghost nodes' mass on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list of masses from jproc:

        jproc = idest->first;

	size_r = idest->second.size();

        vector<float> buf_recv(size_r);

	// cout << "proc " << universe->me << " receives masses from " << jproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_FLOAT, jproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Add the received masses to that of the nodes:

        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            mass[map_ntag[j]] += buf_recv[is];
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

          if (map_ntag[j] != -1) {
            tmp_mass[is] = mass[map_ntag[j]];
          }
        }

	// cout << "proc " << universe->me << " sends mass to " << iproc << endl;
        MPI_Send(tmp_mass.data(), size_s, MPI_FLOAT, iproc, 0, MPI_COMM_WORLD);
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

          if (map_ntag[j] != -1) {
            tmp_mass[is] = mass[map_ntag[j]];
          }
        }

	// cout << "proc " << universe->me << " sends masses to " << jproc << endl;
        MPI_Send(tmp_mass.data(), size_s, MPI_FLOAT, jproc, 0, MPI_COMM_WORLD);
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

	// cout << "proc " << universe->me << " receives size from " << iproc << endl;
        // MPI_Recv(&size_r, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);

	size_r = origin_nshared[iproc].size();
        vector<float> buf_recv(size_r);

	// cout << "proc " << universe->me << " receives masses from " << iproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_FLOAT, iproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Add the received masses to that of the nodes:

        for (int is = 0; is < size_r; is++) {
	  j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            mass[map_ntag[j]] = buf_recv[is];
            // cout << "mass[" << j << "]=" << mass[map_ntag[j]] << "\n";
          }
        }
      }
    }
  }
}

void Grid::reduce_rigid_ghost_nodes() {
  vector<int> tmp_rigid;
  int j, size_r, size_s, jproc;

  // Some ghost nodes' rigid value on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list of rigid from jproc:

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

        vector<int> buf_recv(size_r);

	// cout << "proc " << universe->me << " receives rigids from " << jproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_INT, jproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Apply LOR to the received values of rigid with to that of the nodes:

        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
	    if (buf_recv[is] && !rigid[map_ntag[j]])
	      rigid[map_ntag[j]] = true;
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

        tmp_rigid.assign(size_s, 0); // Reset vector to all 0
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            tmp_rigid[is] = (int) rigid[map_ntag[j]]; // Populate tmp_rigid with the value of rigid
          }
        }

	// cout << "proc " << universe->me << " sends rigids to " << iproc << endl;
	MPI_Send(tmp_rigid.data(),  size_s, MPI_INT, iproc, 0, MPI_COMM_WORLD);
      }
    }
  }

  // Share the updated rigids:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of masses to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

	// cout << "proc " << universe->me << " sends size to " << jproc << endl;
        // MPI_Send(&size_s, 1, MPI_INT, jproc, 0, universe->uworld);

        // Create the list of rigids:

        tmp_rigid.assign(size_s, 0);
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            tmp_rigid[is] = (int) rigid[map_ntag[j]];
          }
        }

	// cout << "proc " << universe->me << " sends rigids to " << jproc << endl;
        MPI_Send(tmp_rigid.data(), size_s, MPI_INT, jproc, 0, MPI_COMM_WORLD);
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

	// cout << "proc " << universe->me << " receives size from " << iproc << endl;
        // MPI_Recv(&size_r, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD,
        //          MPI_STATUS_IGNORE);

	size_r = origin_nshared[iproc].size();
        vector<int> buf_recv(size_r);

	// cout << "proc " << universe->me << " receives rigids from " << iproc << endl;
        MPI_Recv(&buf_recv[0], size_r, MPI_INT, iproc, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // LOR the received rigids to that of the nodes:

        for (int is = 0; is < size_r; is++) {
	  j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
	    if (buf_recv[is] && !rigid[map_ntag[j]])
	      rigid[map_ntag[j]] = true;
          }
        }
      }
    }
  }
}

void Grid::reduce_ghost_nodes(bool reduce_v, bool reduce_forces, bool temp) {
  vector<vector<float>> tmp_vect(universe->nprocs);
  vector<vector<float>> buf_recv_vect(universe->nprocs);
  //vector<float> tmp_mass;
  int j, k, m, size_r, size_s, jproc, nsend;

  nsend = (3 + temp) * (reduce_v + 2*reduce_forces);

  // 1. Pack mass of nodes to be sent back to their owner:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {
      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        jproc = idest->first;
	size_r = idest->second.size();

        buf_recv_vect[jproc].resize(nsend * size_r);
      }
    } else {
      // Check if iproc is listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // me should send the list of masses of all the ghost nodes gotten from
        // iproc to iproc.
        size_s = origin_nshared[iproc].size();

        // Create the list of masses:

        tmp_vect[iproc].resize(nsend * size_s);

	k = 0;

	for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
	    if (reduce_v) {
	      tmp_vect[iproc][k++] = v[m][0];
	      tmp_vect[iproc][k++] = v[m][1];
	      tmp_vect[iproc][k++] = v[m][2];
	    }
            if (reduce_forces) {
              tmp_vect[iproc][k++] = f[m][0];
              tmp_vect[iproc][k++] = f[m][1];
              tmp_vect[iproc][k++] = f[m][2];

              tmp_vect[iproc][k++] = mb[m][0];
              tmp_vect[iproc][k++] = mb[m][1];
              tmp_vect[iproc][k++] = mb[m][2];
            }

	    if (temp) {
	      if (reduce_v)
		tmp_vect[iproc][k++] = T[m];
	      if (reduce_forces) {
		tmp_vect[iproc][k++] = Qint[m];
		tmp_vect[iproc][k++] = Qext[m];
	      }	      
	    }
          } else {
	    error->one(FLERR, "Grid node j does not exist on this CPU.\n");
	  }
        }
      }      
    }
  }

  // 2. New send and receive
  for (int i = 0; i < universe->sendnrecv.size(); i++) {
    if (universe->sendnrecv[i][0] == 0) {
      // Receive
      jproc = universe->sendnrecv[i][1];

      MPI_Recv(&buf_recv_vect[jproc][0], nsend * dest_nshared[jproc].size(), MPI_FLOAT,
               jproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      // Send
      jproc = universe->sendnrecv[i][1];
      MPI_Send(tmp_vect[jproc].data(), nsend * origin_nshared[jproc].size(),
               MPI_FLOAT, jproc, 0, MPI_COMM_WORLD);
    }
  }

  // 3. Add the received masses to that of the nodes:
  for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
       ++idest) {

    k = 0;

    for (int is = 0; is < idest->second.size(); is++) {
      j = idest->second[is];
      if (map_ntag[j] != -1) {
	//mass[map_ntag[j]] += buf_recv_vect[idest->first][is];
	m = map_ntag[j];
	//k = nsend * is;
	if (reduce_v) {
	  v[m][0] += buf_recv_vect[idest->first][k++];
	  v[m][1] += buf_recv_vect[idest->first][k++];
	  v[m][2] += buf_recv_vect[idest->first][k++];
	}
	if (reduce_forces) {
	  f[m][0] += buf_recv_vect[idest->first][k++];
	  f[m][1] += buf_recv_vect[idest->first][k++];
	  f[m][2] += buf_recv_vect[idest->first][k++];

	  mb[m][0] += buf_recv_vect[idest->first][k++];
	  mb[m][1] += buf_recv_vect[idest->first][k++];
	  mb[m][2] += buf_recv_vect[idest->first][k++];
	}

	if (temp) {
	  if (reduce_v)
	    T[m] += buf_recv_vect[idest->first][k++];
	  if (reduce_forces) {
	    Qint[m] += buf_recv_vect[idest->first][k++];
	    Qext[m] += buf_recv_vect[idest->first][k++];
	  }
	}
      } else {
	error->one(FLERR, "Grid node j does not exist on this CPU.\n");
      }
    }
  }

  // Share the updated masses:
  // 1. Pack
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Send the updated list of masses to jproc:

        jproc = idest->first;
        size_s = idest->second.size();

	// cout << "proc " << universe->me << " sends size to " << jproc << endl;

        // Create the list of masses:
        tmp_vect[jproc].resize(nsend * size_s);

	k = 0;
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
	    if (reduce_v) {
	      tmp_vect[jproc][k++] = v[m][0];
	      tmp_vect[jproc][k++] = v[m][1];
	      tmp_vect[jproc][k++] = v[m][2];
	    }
            if (reduce_forces) {
              tmp_vect[jproc][k++] = f[m][0];
              tmp_vect[jproc][k++] = f[m][1];
              tmp_vect[jproc][k++] = f[m][2];

              tmp_vect[jproc][k++] = mb[m][0];
              tmp_vect[jproc][k++] = mb[m][1];
              tmp_vect[jproc][k++] = mb[m][2];
            }

	    if (temp) {
	      if (reduce_v)
		tmp_vect[jproc][k++] = T[m];
	      if (reduce_forces) {
		tmp_vect[jproc][k++] = Qint[m];
		tmp_vect[jproc][k++] = Qext[m];
	      }	      
	    }
          } else {
	    error->one(FLERR, "Grid node j does not exist on this CPU.\n");
	  }
        }
      }
    } else {

      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
        // Receive the updated list of masses from iproc

	// cout << "proc " << universe->me << " receives size from " << iproc << endl;
        buf_recv_vect[iproc].resize(nsend * origin_nshared[iproc].size());
      }
    }
  }


  // 2. New send and receive
  for (int i = 0; i < universe->sendnrecv.size(); i++) {
    if (universe->sendnrecv[i][0] == 0) {
      // Receive
      jproc = universe->sendnrecv[i][1];

      MPI_Send(tmp_vect[jproc].data(), nsend * dest_nshared[jproc].size(),
               MPI_FLOAT, jproc, 0, MPI_COMM_WORLD);
    } else {
      // Send
      jproc = universe->sendnrecv[i][1];
      MPI_Recv(&buf_recv_vect[jproc][0], nsend * origin_nshared[jproc].size(),
               MPI_FLOAT, jproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // 3. Overwrite the nodes' mass with the received values
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc != universe->me) {
      // Check if iproc listed in origin_nshared:

      if (origin_nshared.count(iproc)) {
	k = 0;
        for (int is = 0; is < origin_nshared[iproc].size(); is++) {
	  j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
	    m = map_ntag[j];
            //k = nsend * is;
	    if (reduce_v) {
	      v[m][0] = buf_recv_vect[iproc][k++];
	      v[m][1] = buf_recv_vect[iproc][k++];
	      v[m][2] = buf_recv_vect[iproc][k++];
	    }
            if (reduce_forces) {
              f[m][0] = buf_recv_vect[iproc][k++];
              f[m][1] = buf_recv_vect[iproc][k++];
              f[m][2] = buf_recv_vect[iproc][k++];

              mb[m][0] = buf_recv_vect[iproc][k++];
              mb[m][1] = buf_recv_vect[iproc][k++];
              mb[m][2] = buf_recv_vect[iproc][k++];
            }

            if (temp) {
	      if (reduce_v)
		T[m] = buf_recv_vect[iproc][k++];
              if (reduce_forces) {
                Qint[m] = buf_recv_vect[iproc][k++];
                Qext[m] = buf_recv_vect[iproc][k++];
              }
            }
          } else {
            error->one(FLERR, "Grid node j does not exist on this CPU.\n");
          }
        }
      }
    }
  }
}

void Grid::reduce_ghost_nodes_old(bool only_v, bool temp) {
  vector<float> tmp;
  int j, k, m, size_r, size_s, jproc, nsend;

  if (only_v) {
    nsend = 3 + temp;
  } else {
    nsend = (3 + temp) * 3;
  }

  // Some ghost nodes' mass on the CPU that owns them:
  for (int iproc = 0; iproc < universe->nprocs; iproc++) {
    if (iproc == universe->me) {

      for (auto idest = dest_nshared.cbegin(); idest != dest_nshared.cend();
           ++idest) {
        // Receive the list from jproc:

        jproc = idest->first;

        size_r = idest->second.size();

        vector<float> buf_recv(nsend * size_r);

        // cout << "proc " << universe->me << " receives masses from " << jproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_FLOAT, jproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Add the received data to that of the nodes:

	k = 0;
        for (int is = 0; is < size_r; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
            v[m][0] += buf_recv[k++];
            v[m][1] += buf_recv[k++];
            v[m][2] += buf_recv[k++];
            if (!only_v) {
              f[m][0] += buf_recv[k++];
              f[m][1] += buf_recv[k++];
              f[m][2] += buf_recv[k++];

              mb[m][0] += buf_recv[k++];
              mb[m][1] += buf_recv[k++];
              mb[m][2] += buf_recv[k++];
            }

            if (temp) {
              T[m] += buf_recv[k++];
              if (!only_v) {
                Qint[m] += buf_recv[k++];
                Qext[m] += buf_recv[k++];
              }
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
	k = 0;
        for (int is = 0; is < size_s; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
            tmp[k++] = v[m][0];
            tmp[k++] = v[m][1];
            tmp[k++] = v[m][2];
            if (!only_v) {
              tmp[k++] = f[m][0];
              tmp[k++] = f[m][1];
              tmp[k++] = f[m][2];

              tmp[k++] = mb[m][0];
              tmp[k++] = mb[m][1];
              tmp[k++] = mb[m][2];
            }

	    if (temp) {
	      tmp[k++] = T[m];
	      if (!only_v) {
		tmp[k++] = Qint[m];
		tmp[k++] = Qext[m];
	      }	      
	    }
          }
        }

        // cout << "proc " << universe->me << " sends mass to " << iproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_FLOAT, iproc, 0,
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
	k = 0;
        for (int is = 0; is < size_s; is++) {
          j = idest->second[is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
            tmp[k++] = v[m][0];
            tmp[k++] = v[m][1];
            tmp[k++] = v[m][2];
            if (!only_v) {
              tmp[k++] = f[m][0];
              tmp[k++] = f[m][1];
              tmp[k++] = f[m][2];

              tmp[k++] = mb[m][0];
              tmp[k++] = mb[m][1];
              tmp[k++] = mb[m][2];
            }

	    if (temp) {
	      tmp[k++] = T[m];
	      if (!only_v) {
		tmp[k++] = Qint[m];
		tmp[k++] = Qext[m];
	      }	      
	    }
          }
        }

        // cout << "proc " << universe->me << " sends masses to " << jproc <<
        // endl;
        MPI_Send(tmp.data(), nsend * size_s, MPI_FLOAT, jproc, 0,
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
        vector<float> buf_recv(nsend * size_r);

        // cout << "proc " << universe->me << " receives masses from " << iproc
        // << endl;
        MPI_Recv(&buf_recv[0], nsend * size_r, MPI_FLOAT, iproc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Update the local data with the received values
	k = 0;
        for (int is = 0; is < size_r; is++) {
          j = origin_nshared[iproc][is];

          if (map_ntag[j] != -1) {
            m = map_ntag[j];
            //k = nsend * is;
            v[m][0] = buf_recv[k++];
            v[m][1] = buf_recv[k++];
            v[m][2] = buf_recv[k++];
            if (!only_v) {
              f[m][0] = buf_recv[k++];
              f[m][1] = buf_recv[k++];
              f[m][2] = buf_recv[k++];

              mb[m][0] = buf_recv[k++];
              mb[m][1] = buf_recv[k++];
              mb[m][2] = buf_recv[k++];
            }

            if (temp) {
              T[m] = buf_recv[k++];
              if (!only_v) {
                Qint[m] = buf_recv[k++];
                Qext[m] = buf_recv[k++];
              }
            }
          }
        }
      }
    }
  }
}
