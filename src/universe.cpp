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

#include <mpi.h>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <universe.h>
#include <version.h>
#include <domain.h>
#include <error.h>

using namespace std;

static vector<int> tile2d(int); ///< Determines the tiling of p procs in 2D
static vector<int> tile3d(int); ///< Determines the tiling of p proces in 3D

struct boundsize { 
    double dl;
    int rank;
};

Universe::Universe(MPM *mpm, MPI_Comm communicator) : Pointers(mpm)
{
  uworld = communicator;
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);
}

Universe::~Universe()
{
  // if (uworld != uorig) MPI_Comm_free(&uworld);
}

void Universe::set_proc_grid() {
  int dim = domain->dimension;

  if (dim != 1 && dim !=2 && dim !=3) {
    error->all(FLERR, "Error in Universe::set_proc_grid(): invalid dimension: " + to_string(dim) + ".\n");
  }

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  procgrid[0] = 1;
  procgrid[1] = 1;
  procgrid[2] = 1;

  myloc[0] = 0;
  myloc[1] = 0;
  myloc[2] = 0;

  if (nprocs > 1 && dim == 1) {
    // Easy bit, all procs are in line:
    procgrid[0] = nprocs;

    myloc[0] = me;
  }

  if (nprocs > 1 && dim == 2) {
    // Determine the smallest dimension:
    double l[2] = {domain->boxhi[0] - domain->boxlo[0],
		   domain->boxhi[1] - domain->boxlo[1]};

    if (l[0]<1.0e-10 || l[1]<1.0e-10) {
      error->all(FLERR, "Error: the domain has a size in at least one direction that is 0: Lx=" + to_string(l[0]) + ", Ly=" + to_string(l[1]) + ".\n");
    }

    // length >= width
    int length = 1; // Which direction represents the domain's length
    int width = 0;  // Which direction represents the domain's width

    if (l[0] > l[1]) {
      width = 1;
      length = 0;
    }

    vector<int> tile = tile2d(nprocs);
    procgrid[width] = tile[0];
    procgrid[length] = tile[1];

    myloc[0] = me % procgrid[0];
    myloc[1] = me / procgrid[0];
  }

  if (nprocs > 1 && dim == 3) {
    // Determine the smallest dimension:
    vector<boundsize> l = {{domain->boxhi[0] - domain->boxlo[0], 0},
			   {domain->boxhi[1] - domain->boxlo[1], 1},
			   {domain->boxhi[2] - domain->boxlo[2], 2}};

    if (l[0].dl<1.0e-10 || l[1].dl<1.0e-10 || l[2].dl<1.0e-10) {
      error->all(FLERR, "Error: the domain has a size in at least one direction that is 0: Lx=" + to_string(l[0].dl) + ", Ly=" + to_string(l[1].dl) + ", Lz=" + to_string(l[2].dl) + ".\n");
    }

    // Sort values in l according to their dl:
    sort(l.begin(), l.end(),
	 [](boundsize const &a, boundsize const &b) { return a.dl < b.dl; });

    vector<int> tile = tile3d(nprocs);
    procgrid[l[0].rank] = tile[0];
    procgrid[l[1].rank] = tile[1];
    procgrid[l[2].rank] = tile[2];
    
    myloc[2] = me / (procgrid[0]*procgrid[1]);
    myloc[1] = me / procgrid[0] - myloc[2]*procgrid[1];
    myloc[0] = me - (myloc[1] + myloc[2]*procgrid[1])*procgrid[0];
  }


  domain->set_local_box();


  procneigh[0][0] = procneigh[0][1] = 0;
  procneigh[1][0] = procneigh[1][1] = 0;
  procneigh[2][0] = procneigh[2][1] = 0;

  if (dim == 1) {
    if (sublo[0] > boxlo[0] + 1.0e-12) procneigh[0][0] = me - 1;
    else procneigh[0][0] = -1;
    if (subhi[0] < boxhi[0] - 1.0e-12) procneigh[0][1] = me + 1;
    else procneigh[0][1] = -1;
  }

  if (dim == 2) {
    if (sublo[0] > boxlo[0] + 1.0e-12)
      procneigh[0][0] = me - 1;
    else
      procneigh[0][0] = -1;
    if (subhi[0] < boxhi[0] - 1.0e-12)
      procneigh[0][1] = me + 1;
    else
      procneigh[0][1] = -1;
    if (sublo[1] > boxlo[1] + 1.0e-12)
      procneigh[1][0] = me - procgrid[0];
    else
      procneigh[1][0] = -1;
    if (subhi[1] < boxhi[1] - 1.0e-12)
      procneigh[1][1] = me + procgrid[0];
    else
      procneigh[1][1] = -1;
  }

  if (dim == 3) {
    if (sublo[0] > boxlo[0] + 1.0e-12)
      procneigh[0][0] = me - 1;
    else
      procneigh[0][0] = -1;
    if (subhi[0] < boxhi[0] - 1.0e-12)
      procneigh[0][1] = me + 1;
    else
      procneigh[0][1] = -1;
    if (sublo[1] > boxlo[1] + 1.0e-12)
      procneigh[1][0] = me - procgrid[0];
    else
      procneigh[1][0] = -1;
    if (subhi[1] < boxhi[1] - 1.0e-12)
      procneigh[1][1] = me + procgrid[0];
    else
      procneigh[1][1] = -1;
    if (sublo[2] > boxlo[2] + 1.0e-12)
      procneigh[2][0] = me - procgrid[0]*procgrid[1];
    else
      procneigh[2][0] = -1;
    if (subhi[2] < boxhi[2] - 1.0e-12)
      procneigh[2][1] = me + procgrid[0]*procgrid[1];
    else
      procneigh[2][1] = -1;
  }

#ifdef DEBUG
    cout << "proc " << universe->me << "\tprocneigh = [[" << procneigh[0][0] << "," << procneigh[0][1] <<"],[" << procneigh[1][0] << "," << procneigh[1][1] <<"],[" << procneigh[2][0] << "," << procneigh[2][1] << "].\n";
#endif

    // Set send and receive pattern:
    if (domain->dimension >= 2) {
      int jproc;
      // Step 1:
      if (me % 2 == 1) {
	jproc = procneigh[0][0];
	if (jproc != -1) {
	  sendnrecv.push_back({1, jproc, 1}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 1}); // me receives from jproc
	}
      }

      if (me % 2 == 0) {
	jproc = procneigh[0][1];
	if (jproc != -1) {
	  sendnrecv.push_back({0, jproc, 1}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 1}); // me sends to jproc
	}
      }

      // Step 2:
      if (me % 2 == 0) {
	jproc = procneigh[0][0];
	if (jproc != -1) {
	  sendnrecv.push_back({1, jproc, 2}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 2}); // me receives from jproc
	}
      }

      if (me % 2 == 1) {
	jproc = procneigh[0][1];
	if (jproc != -1) {
	  sendnrecv.push_back({0, jproc, 2}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 2}); // me sends to jproc
	}
      }


      // Step 3:
      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][0];
	if (jproc != -1) {
	  sendnrecv.push_back({1, jproc, 3}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 3}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][1];
	if (jproc != -1) {
	  sendnrecv.push_back({0, jproc, 3}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 3}); // me sends to jproc
	}
      }


      // Step 4:
      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][0];
	if (jproc != -1) {
	  sendnrecv.push_back({1, jproc, 4}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 4}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][1];
	if (jproc != -1) {
	  sendnrecv.push_back({0, jproc, 4}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 4}); // me sends to jproc
	}
      }


      // Step 5:
      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][0] + 1;
	if (jproc > 0 && jproc / procgrid[0] < me / procgrid[0]) {
	  sendnrecv.push_back({1, jproc, 5}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 5}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][1] - 1;
	if (jproc > -1 && me / procgrid[0] < jproc / procgrid[0]) {
	  sendnrecv.push_back({0, jproc, 5}); // me receives from jproc
	  sendnrecv.push_back({1, jproc,5}); // me sends to jproc
	}
      }

      // Step 6:
      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][0] + 1;
	if (jproc > 0 && jproc / procgrid[0] < me / procgrid[0] && procneigh[0][1] != -1) {
	  sendnrecv.push_back({1, jproc, 6}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 6}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][1] - 1;
	if (jproc != -1 && me / procgrid[0] < jproc / procgrid[0]) {
	  sendnrecv.push_back({0, jproc, 6}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 6}); // me sends to jproc
	}
      }


      // Step 7:
      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][0] - 1;
	if (jproc > -1 && me / procgrid[0] - jproc / procgrid[0] == 1) {
	  sendnrecv.push_back({1, jproc, 7}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 7}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][1] + 1;
	if (jproc > 0 && jproc / procgrid[0] - me / procgrid[0] == 1) {
	  sendnrecv.push_back({0, jproc, 7}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 7}); // me sends to jproc
	}
      }


      // Step 8:
      if ((me / procgrid[0]) % 2 == 0) {
	jproc = procneigh[1][0] - 1;
	if (jproc > -1 && me / procgrid[0] - jproc / procgrid[0] == 1) {
	  sendnrecv.push_back({1, jproc, 8}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 8}); // me receives from jproc
	}
      }

      if ((me / procgrid[0]) % 2 == 1) {
	jproc = procneigh[1][1] + 1;
	if (jproc > 0 && jproc / procgrid[0] - me / procgrid[0] == 1) {
	  sendnrecv.push_back({0, jproc, 8}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 8}); // me sends to jproc
	}
      }

      if (me == 0)
	cout << "End set communication pattern\n";
    }

    if (domain->dimension == 3) {
      int jproc;
      // Step 9:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0];
	if (jproc > -1) {
	  sendnrecv.push_back({1, jproc, 9}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 9}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1];
	if (jproc > -1) {
	  sendnrecv.push_back({0, jproc, 9}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 9}); // me sends to jproc
	}
      }

      // Step 10:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0];
	if (jproc > -1) {
	  sendnrecv.push_back({1, jproc, 10}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 10}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1];
	if (jproc > -1) {
	  sendnrecv.push_back({0, jproc, 10}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 10}); // me sends to jproc
	}
      }

      // Step 11:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] - procgrid[0];
	if (jproc > -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({1, jproc, 11}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 11}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] + procgrid[0];
	if (jproc > -1 && procneigh[1][1] != -1 && procneigh[2][1] != -1) {
	  sendnrecv.push_back({0, jproc, 11}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 11}); // me sends to jproc
	}
      }

      // Step 12:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] - procgrid[0];
	if (jproc > -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({1, jproc, 12}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 12}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] + procgrid[0];
	if (jproc > -1 && procneigh[1][1] != -1 && procneigh[2][1] != -1) {
	  sendnrecv.push_back({0, jproc, 12}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 12}); // me sends to jproc
	}
      }

      // Step 13:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] - 1;
	if (jproc > -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 13}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 13}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] + 1;
	if (jproc > -1 && procneigh[2][1] != -1 && procneigh[0][1] != -1) {
	  sendnrecv.push_back({0, jproc, 13}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 13}); // me sends to jproc
	}
      }

      // Step 14:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] - 1;
	if (jproc > -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 14}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 14}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] + 1;
	if (jproc > -1 && procneigh[0][1] != -1 && procneigh[2][1] != -1) {
	  sendnrecv.push_back({0, jproc, 14}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 14}); // me sends to jproc
	}
      }

      // Step 15:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] - procgrid[0] - 1;
	if (jproc > -1 && procneigh[1][0] != -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 15}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 15}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] + procgrid[0] + 1;
	if (jproc > procgrid[0] && procneigh[1][1] != -1 && procneigh[0][1] != -1) {
	  sendnrecv.push_back({0, jproc, 15}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 15}); // me sends to jproc
	}
      }

      // Step 16:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] - procgrid[0] - 1;
	if (jproc > -1 && procneigh[1][0] != -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 16}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 16}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] + procgrid[0] + 1;
	if (jproc > procgrid[0] && procneigh[1][1] != -1 && procneigh[0][1] != -1) {
	  sendnrecv.push_back({0, jproc, 16}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 16}); // me sends to jproc
	}
      }

      // Step 17:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] + procgrid[0] - 1;
	if (jproc > -1 && procneigh[1][1] != -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 17}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 17}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] - procgrid[0] + 1;
	if (jproc > -1 && procneigh[1][0] != -1 && procneigh[0][1] != -1) {
	  sendnrecv.push_back({0, jproc, 18}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 17}); // me sends to jproc
	}
      }

      // Step 18:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] + procgrid[0] - 1;
	if (jproc > procgrid[0] - 2 && procneigh[1][1] != -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({1, jproc, 18}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 18}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] - procgrid[0] + 1;
	if (jproc > -1 && procneigh[1][0] != -1 && procneigh[0][1] != -1) {
	  sendnrecv.push_back({0, jproc, 18}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 18}); // me sends to jproc
	}
      }

      // Step 19:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] + procgrid[0];
	if (jproc > -1 && procneigh[1][1] != -1) {
	  sendnrecv.push_back({1, jproc, 19}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 19}); // me receives from jproc
	}
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] - procgrid[0];
	if (jproc > -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({0, jproc, 19}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 19}); // me sends to jproc
	}
      }

      // Step 20:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] + procgrid[0];
        if (jproc > procgrid[0] - 1 && procneigh[1][1] != -1) {
          sendnrecv.push_back({1, jproc, 20}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 20}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] - procgrid[0];
	if (jproc > -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({0, jproc, 20}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 20}); // me sends to jproc
	}
      }

      // Step 21:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] + procgrid[0] + 1;
	if (jproc > -1 && procneigh[0][1] != -1 && procneigh[1][1] != -1) {
          sendnrecv.push_back({1, jproc, 21}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 21}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] - procgrid[0] - 1;
	if (jproc > -1 && procneigh[0][0] != -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({0, jproc, 21}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 21}); // me sends to jproc
	}
      }

      // Step 22:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] + procgrid[0] + 1;
	if (jproc > procgrid[0] && procneigh[0][1] != -1 && procneigh[1][1] != -1) {
          sendnrecv.push_back({1, jproc, 22}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 22}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] - procgrid[0] - 1;
	if (jproc > -1 && procneigh[0][0] != -1 && procneigh[1][0] != -1) {
	  sendnrecv.push_back({0, jproc, 22}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 22}); // me sends to jproc
	}
      }

      // Step 23:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] + 1;
	if (jproc > -1 && procneigh[0][1] != -1) {
          sendnrecv.push_back({1, jproc, 23}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 23}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] - 1;
	if (jproc > -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({0, jproc, 23}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 23}); // me sends to jproc
	}
      }

      // Step 24:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] + 1;
	if (jproc > 0 && procneigh[0][1] != -1) {
          sendnrecv.push_back({1, jproc, 24}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 24}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] - 1;
	if (jproc > -1 && procneigh[0][0] != -1) {
	  sendnrecv.push_back({0, jproc, 24}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 24}); // me sends to jproc
	}
      }

      // Step 25:
      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][0] - procgrid[0] + 1;
	if (jproc > -1 && procneigh[0][1] != -1 && procneigh[1][0] != -1) {
          sendnrecv.push_back({1, jproc, 25}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 25}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][1] + procgrid[0] - 1;
	if (jproc > -1 && procneigh[0][0] != -1 && procneigh[1][1] != -1 && procneigh[2][1] != -1) {
	  sendnrecv.push_back({0, jproc, 25}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 25}); // me sends to jproc
	}
      }

      // Step 26:
      if (myloc[2] % 2 == 0) {
	jproc = procneigh[2][0] - procgrid[0] + 1;
	if (jproc > -1 && procneigh[0][1] != -1 && procneigh[1][0] != -1) {
          sendnrecv.push_back({1, jproc, 26}); // me sends to jproc
	  sendnrecv.push_back({0, jproc, 26}); // me receives from jproc
        }
      }

      if (myloc[2] % 2 == 1) {
	jproc = procneigh[2][1] + procgrid[0] - 1;
	if (jproc > procgrid[0] - 2 && procneigh[0][0] != -1 && procneigh[1][1] != -1) {
	  sendnrecv.push_back({0, jproc, 26}); // me receives from jproc
	  sendnrecv.push_back({1, jproc, 26}); // me sends to jproc
	}
      }
    }

    if (domain->dimension == 1) {
      error->all(FLERR, "New partitioning not supported for dimensions 1 yet!\n");
    }

    // For debug purposes only:
    for (int p=0; p<nprocs; p++) {
      if (p==me) {
	cout << "proc " << me << " sendnrecv=[";
	for (int i=0; i<sendnrecv.size(); i++) {
	  cout << "[" << sendnrecv[i][0] << ", " << sendnrecv[i][1] << ", " << sendnrecv[i][2] << "],";
	}
	cout << "]\tmyloc=[" << myloc[0] << "," << myloc[1] << "," << myloc[2] << "]\n";
      }
      MPI_Barrier(uworld);
    }
}

vector<int> tile2d(int p) {
  vector<int> result;
  int n = (int) floor(sqrt(p));

  while ((p/n)*n != p) {
    n -= 1;
  }

  result.push_back(n);
  result.push_back(p/n);

  return result;
}


vector<int> tile3d(int p) {
  vector<int> result;

  int n = round(pow(p,1.0/3.0));

  while ((p/n)*n != p) {
    n -= 1;
  }

  result.push_back(n);
  vector<int> ml = tile2d(p/n);

  result.push_back(ml[0]);
  result.push_back(ml[1]);
  sort(result.begin(), result.end(), greater<int>());

  return result;
}
