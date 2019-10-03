#include <iostream>
#include "group.h"
#include "delete_particles.h"
#include "domain.h"
#include "error.h"

using namespace std;

DeleteParticles::DeleteParticles(MPM *mpm) : Pointers(mpm) {}

Var DeleteParticles::command(vector<string> args) {
  cout << "In DeleteParticles::command()" << endl;

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
  // for (int i=0; i<ns; i++){
  //   if (dlist[i]!=NULL) {
  //     int np = domain->solids[i]->np;
  //     double vtot = domain->solids[i]->vtot;

  //     int k = 0;
  //     while (k < np) {
  // 	if (dlist[i][k]) {
  // 	  vtot -= domain->solids[i]->vol[k];
  // 	  domain->solids[i]->copy_particle(np-1,k);
  // 	  dlist[i][k] = dlist[i][np-1];
  // 	  np--;
  // 	} else k++;
  //     }
  //     domain->solids[i]->vtot = vtot;
  //     cout << "Deleting " << domain->solids[i]->np - np << " particles from solid " << domain->solids[i]->id << endl;
  //     cout << "Solid " << domain->solids[i]->id << " new total volume = " << domain->solids[i]->vtot << endl;
  //     domain->solids[i]->np = np;
  //   }
  // }
  
  return Var(0);
}

void DeleteParticles::delete_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0) error->all(FLERR, "Error: region " + args[2] + " unknown.\n");

  int ns = domain->solids.size();
  int np;
  for(int is=0; is<ns; is++){
    if ((isolid < 0) || (is == isolid)) {
      np = domain->solids[isolid]->np;
      dlist[is] = new int[np];

      for(int ip=0; ip<np; ip++){
	if (domain->regions[iregion]->inside(domain->solids[is]->x0[ip][0],
					     domain->solids[is]->x0[ip][1],
					     domain->solids[is]->x0[ip][2])==1) dlist[is][ip] = 1;
	else dlist[is][ip] = 0;
      }  
    }
  }
}
