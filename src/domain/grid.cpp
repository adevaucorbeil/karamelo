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
  nsolids = 0;
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
  int nsolids = this->nsolids;

  Kokkos::View<tagint*> ntag = this->ntag;
  Kokkos::View<int*> nowner = this->nowner;

  Kokkos::View<Vector3d*> x = this->x;
  Kokkos::View<Vector3d*> x0 = this->x0;

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

    ntag[l] = nz_global*ny_global*(i+noffsetlo_0) + nz_global*(j+noffsetlo_1) + k+noffsetlo_2;
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
  }
}

void Grid::setup(string cs){
  cellsize = input->parsev(cs);
  if (universe->me == 0) {
    cout << "Set grid cellsize to " << cellsize << endl;
  }

  // Check if the size of the non-shared views (mass, v, v_update, ...)
  // is appropriate for the number of solids present:

  if (update->method->slip_contacts && mass.extent(0) < nsolids) {
    int nn = mass.extent(1);
    Kokkos::resize(v,        nsolids, nn);
    Kokkos::resize(v_update, nsolids, nn);
    Kokkos::resize(mb,       nsolids, nn);
    Kokkos::resize(f,        nsolids, nn);
    Kokkos::resize(mass,     nsolids, nn);
    if (update->method->anti_volumetric_locking)
      Kokkos::resize(vol,    nsolids, nn);

    Kokkos::resize(rigid,    nsolids, nn);
    Kokkos::parallel_for(__PRETTY_FUNCTION__,
     Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nsolids, nn}),
    KOKKOS_LAMBDA (const int &is, const int &in)
    {
      rigid(is, in) = false;
    });

    Kokkos::resize(T,        nsolids, nn);
    Kokkos::resize(T_update, nsolids, nn);
    Kokkos::resize(Qext,     nsolids, nn);
    Kokkos::resize(Qint,     nsolids, nn);
    Kokkos::resize(normal,   nsolids, nn);
  }
}

void Grid::grow(int nn){
  //nnodes_local = nn;

  ntag   = Kokkos::View<tagint*>("ntag",   nn);
  nowner = Kokkos::View<int*>   ("nowner", nn);

  x0       = Kokkos::View<Vector3d*>("x0",       nn);
  x        = Kokkos::View<Vector3d*>("x",        nn);

  int ns = nsolids >= 1 && update->method->slip_contacts? nsolids: 1;

  v        = Kokkos::View<Vector3d**>("v",        ns, nn);
  v_update = Kokkos::View<Vector3d**>("v_update", ns, nn);
  mb       = Kokkos::View<Vector3d**>("mb",       ns, nn);
  f        = Kokkos::View<Vector3d**>("f",        ns, nn);

  mass = Kokkos::View<float**>("mass", ns, nn);
  if (update->method->anti_volumetric_locking)
    vol = Kokkos::View<float**>("vol", ns, nn);
  mask = Kokkos::View<int*>   ("mask", nn);

  Kokkos::View<int*> mask = this->mask;

  Kokkos::parallel_for(__PRETTY_FUNCTION__, nn,
  KOKKOS_LAMBDA(int i)
  {
    mask[i] = 1;
  });

  ntype = Kokkos::View<Vector3i*>("ntype", nn);
  rigid = Kokkos::View<bool**>   ("rigid", ns, nn);
  Kokkos::parallel_for(__PRETTY_FUNCTION__,
     Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ns, nn}),
  KOKKOS_LAMBDA (const int &is, const int &in)
  {
    rigid(is, in) = false;
  });

  T        = Kokkos::View<float**>("T",        ns, nn);
  T_update = Kokkos::View<float**>("T_update", ns, nn);
  Qext     = Kokkos::View<float**>("Qext",     ns, nn);
  Qint     = Kokkos::View<float**>("Qint",     ns, nn);

  normal   = Kokkos::View<Vector3d**>("normal", ns, nn);
}
