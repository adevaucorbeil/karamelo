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

#include "delete_particles.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "universe.h"
#include "var.h"
#include <iostream>

using namespace std;

DeleteParticles::DeleteParticles(MPM *mpm) : Pointers(mpm) {}

Var DeleteParticles::command(vector<string> args) {
  // cout << "In DeleteParticles::command()" << endl;

  if (args.size() < 3) error->all(FLERR, "Illegal delete_particles command\n");

  int ns = domain->solids.size();
  dlist = new int*[ns];

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) error->all(FLERR, "Error: solid " + args[0] + " unknown.\n");
  }

  if (args[1].compare("region")==0) delete_region(args, isolid);
  else error->all(FLERR, "Error: use of illegal keyword for delete_particles command: " + args[1] + "\n");

  // Delete particles flagged in dlist:
  for (int i = 0; i < ns; i++) {
    if (dlist[i] != nullptr) {
      int np_local = domain->solids[i]->np_local;
      double vtot_local = 0;
      double mtot_local = 0;
      for (int ip=0; ip<np_local; ip++)
	{
	  vtot_local += domain->solids[i]->vol[ip];
	  mtot_local += domain->solids[i]->mass[ip];
	}

      int k = 0;
      while (k < np_local) {
        if (dlist[i][k]) {
          vtot_local -= domain->solids[i]->vol[k];
          mtot_local -= domain->solids[i]->mass[k];
          domain->solids[i]->copy_particle(np_local - 1, k);
          dlist[i][k] = dlist[i][np_local - 1];
          np_local--;
        } else
          k++;
      }

      int np_local_reduced = 0;
      MPI_Allreduce(&np_local, &np_local_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);
      MPI_Allreduce(&vtot_local, &domain->solids[i]->vtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
      MPI_Allreduce(&mtot_local, &domain->solids[i]->mtot, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

      if (universe->me == 0) {
	cout << "Deleting " << domain->np_total - np_local_reduced
	     << " particles from solid " << domain->solids[i]->id << endl;
	cout << "Solid " << domain->solids[i]->id
	     << " new total volume = " << domain->solids[i]->vtot << endl;
      }

      domain->solids[i]->np_local = np_local;
      domain->np_total = np_local_reduced;

    }
  }

  return Var(0);
}

void DeleteParticles::delete_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0)
    error->all(FLERR, "Error: region " + args[2] + " unknown.\n");

  int ns = domain->solids.size();
  int np_local;
  for (int is = 0; is < ns; is++) {
    if ((isolid < 0) || (is == isolid)) {
      np_local = domain->solids[is]->np_local;
      dlist[is] = new int[np_local];

      for (int ip = 0; ip < np_local; ip++) {
        if (domain->regions[iregion]->inside(
                domain->solids[is]->x0[ip][0], domain->solids[is]->x0[ip][1],
                domain->solids[is]->x0[ip][2]) == 1)
          dlist[is][ip] = 1;
        else
          dlist[is][ip] = 0;
      }
    }
  }
}
