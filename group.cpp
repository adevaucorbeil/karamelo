#include "group.h"
#include "domain.h"
#include "region.h"
#include "solid.h"
#include <Eigen/Eigen>

#define MAX_GROUP 32

using namespace std;
using namespace Eigen;

Group::Group(MPM *mpm) : Pointers(mpm)
{
  names = new string[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++) names[i] = "";
  for (int i = 0; i < MAX_GROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++) inversemask[i] = bitmask[i] ^ ~0;

  // create "all" group
  names[0] = "all";
  ngroup = 1;

}

Group::~Group()
{
  delete [] names;
  delete [] bitmask;
  delete [] inversemask;
}

/* ----------------------------------------------------------------------
   assign atoms to a new or existing group
------------------------------------------------------------------------- */

void Group::assign(vector<string> args)
{
  if (args.size() < 4) {
    cout << "Error: too few arguments for group: requires at least 4 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  // style = region
  // add to group if atom is in region

  if (args[1].compare("region") == 0) {
    // Look for the region ID (if exists):
    int iregion = domain->find_region(args[2]);
    if (iregion == -1) {
      cout << "Error: could not find region " << args[2] << endl;
      exit(1);
    }

    // find group in existing list
    // add a new group if igroup = -1

    int igroup = find(args[0]);

    if (igroup == -1) {
      if (ngroup == MAX_GROUP) {
	cout << "Too many groups" << endl;
      }
      igroup = find_unused();
      names[igroup] = args[0];
      ngroup++;
    }

    int bit = bitmask[igroup];

    if (args[3].compare("all") == 0) {

      /* For all particles of all solids, check if they are in the region.
	 If so asign them the right mask */
      for (int isolid = 0; isolid < domain->solids.size(); isolid++) {

	Eigen::Vector3d *x = domain->solids[isolid]->x;
	int *mask = domain->solids[isolid]->mask;
	int n = 0;

	for (int ip = 0; ip < domain->solids[isolid]->np; ip++) {
	  if (domain->regions[iregion]->match(x[ip][0],x[ip][1],x[ip][2])) {
	    mask[ip] |= bit;
	    n++;
	  }
	}
	cout << n << " particles from solid " << domain->solids[isolid]->id << " found" << endl;
      }
    } else if (args[3].compare("solid") == 0) {

      for (int i=4; i<args.size(); i++) {
	int isolid = domain->find_solid(args[i]);
	if (isolid == -1) {
	  cout << "Error: cannot find solid with ID " << args[i] << endl;
	  exit(1);
	}

	Eigen::Vector3d *x = domain->solids[isolid]->x;
	int *mask = domain->solids[isolid]->mask;
	int n = 0;

	for (int ip = 0; ip < domain->solids[isolid]->np; ip++) {
	  if (domain->regions[iregion]->match(x[ip][0],x[ip][1],x[ip][2]))
	    mask[ip] |= bit;
	    n++;
	}
	cout << n << " particles from solid " << domain->solids[isolid]->id << " found" << endl;
      }

    } else {
      cout << "Error: unknown keyword in group command: " << args[3] << endl;
      exit(1);
    }
  }
}

/* ----------------------------------------------------------------------
   return group index if name matches existing group, -1 if no such group
------------------------------------------------------------------------- */

int Group::find(string name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (name.compare(names[igroup]) == 0) return igroup;
  return -1;
}

/* ----------------------------------------------------------------------
   return index of first available group
   should never be called when group limit has been reached
------------------------------------------------------------------------- */

int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == "") return igroup;
  return -1;
}
