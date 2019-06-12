#include <iostream>
#include "group.h"
#include "delete_particles.h"
#include "domain.h"

using namespace std;

DeleteParticles::DeleteParticles(MPM *mpm) : Pointers(mpm) {}

Var DeleteParticles::command(vector<string> args) {
  cout << "In DeleteParticles::command()" << endl;

  if (args.size() < 3) {
    cout << "Illegal delete_particles command" << endl;
    exit(1);
  }

  int ns = domain->solids.size();
  dlist = new int*[ns];

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) {
      cout << "Error: solid " << args[0] << " unknown.\n";
      exit(1);
    }
  }

  if (args[1].compare("region")==0) delete_region(args, isolid);
  else {
    cout << "Error: use of illegal keyword for delete_particles command: " << args[1] << endl;
    exit(1);
  }


  // Delete particles flagged in dlist:
  for (int i=0; i<ns; i++){
    if (dlist[i]!=NULL) {
      int np = domain->solids[i]->np;

      int k = 0;
      while (k < np) {
	if (dlist[i][k]) {
	  domain->solids[i]->copy_particle(np-1,k);
	  dlist[i][k] = dlist[i][np-1];
	  np--;
	} else k++;
      }
      cout << "Deleting " << domain->solids[i]->np - np << " particles from solid " << domain->solids[i]->id << endl;
      domain->solids[i]->np = np;
    }
  }
  
  return Var(0);
}

void DeleteParticles::delete_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0) {
    cout << "Error: region " << args[2] << " unknown.\n";
    exit(1);
  }

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
