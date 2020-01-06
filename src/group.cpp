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

#include <Eigen/Eigen>
#include <string>
#include "group.h"
#include "domain.h"
#include "region.h"
#include "solid.h"
#include "error.h"
#include "universe.h"

#define MAX_GROUP 32

using namespace std;
using namespace Eigen;

Group::Group(MPM *mpm) : Pointers(mpm)
{
  names       = new string[MAX_GROUP];
  bitmask     = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  pon         = new string[MAX_GROUP];
  solid       = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++)
    names[i] = "";
  for (int i = 0; i < MAX_GROUP; i++)
    bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++)
    inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++)
    pon[i] = "all";
  for (int i = 0; i < MAX_GROUP; i++)
    solid[i] = -1;

  // create "all" group
  names[0] = "all";
  ngroup   = 1;
}

Group::~Group()
{
  delete[] names;
  delete[] bitmask;
  delete[] inversemask;
  delete[] pon;
  delete[] solid;
}

/* ----------------------------------------------------------------------
   assign atoms to a new or existing group
------------------------------------------------------------------------- */

void Group::assign(vector<string> args)
{
  if (args.size() < 4) {
    error->all(FLERR,"Error: too few arguments for group: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  // find group in existing list
  // add a new group if igroup = -1

  int igroup = find(args[0]);

  if (igroup == -1)
  {
    if (ngroup == MAX_GROUP)
    {
      cout << "Too many groups" << endl;
    }
    igroup        = find_unused();
    names[igroup] = args[0];
    ngroup++;
  }

  int bit = bitmask[igroup];

  // Operates on particles or nodes:
  if (args[1].compare("particles") == 0)
  {
    pon[igroup] = "particles";
  }
  else if (args[1].compare("nodes") == 0)
  {
    pon[igroup] = "nodes";
  } else {
    error->all(FLERR,"Error: do not understand keyword " + args[1] + ", \"particles\" or \"nodes\" expected.\n");
  }

  // style = region
  // add to group if atom is in region

  if (args[2].compare("region") == 0)
  {
    // Look for the region ID (if exists):
    int iregion = domain->find_region(args[3]);
    if (iregion == -1)
      {
      error->all(FLERR, "Error: could not find region " + args[3] + ".\n");
      }

    if (args[4].compare("all") == 0)
      {

	/* For all particles of all solids, check if they are in the region.
	   If so asign them the right mask */
	solid[igroup] = -1; // Consider all solids

	for (int isolid = 0; isolid < domain->solids.size(); isolid++)
	  {

	    vector<Eigen::Vector3d> *x;
	    int nmax;
	    vector<int> *mask;

	    if (pon[igroup].compare("particles") == 0)
	      {
		x    = &domain->solids[isolid]->x;
		nmax = domain->solids[isolid]->np_local;
		mask = &domain->solids[isolid]->mask;
		cout << "Solid has " << domain->solids[isolid]->np << " particles"
		     << endl;
	      }
	    else
	      {
		x    = &domain->solids[isolid]->grid->x;
		nmax = domain->solids[isolid]->grid->nnodes_local
		  + domain->solids[isolid]->grid->nnodes_ghost;
		mask = &domain->solids[isolid]->grid->mask;
		cout << "Grid has " << domain->solids[isolid]->grid->nnodes
		     << " nodes" << endl;
	      }

	    int n = 0;

	    for (int ip = 0; ip < nmax; ip++)
	      {
		if (domain->regions[iregion]->match((*x)[ip][0],(*x)[ip][1],(*x)[ip][2]))
		  {
		    (*mask)[ip] |= bit;
		    n++;
		  }
	      }

	    int n_tot = 0;
	    MPI_Allreduce(&n,&n_tot,1,MPI_INT,MPI_SUM,universe->uworld);
	
	    cout << n_tot << " " << pon[igroup] << " from solid "
		 << domain->solids[isolid]->id << " found" << endl;
	  }
      }
    else if (args[4].compare("solid") == 0)
      {

	for (int i=5; i<args.size(); i++)
	  {
	    solid[igroup] = domain->find_solid(args[i]);
	    if (solid[igroup] == -1) {
	      error->all(FLERR, "Error: cannot find solid with ID " + args[i] + ".\n");
	    }

	    vector<Eigen::Vector3d> *x;
	    int nmax;
	    vector<int> *mask;

	    if (pon[igroup].compare("particles") == 0)
	      {
		x    = &domain->solids[solid[igroup]]->x;
		nmax = domain->solids[solid[igroup]]->np_local;
		mask = &domain->solids[solid[igroup]]->mask;
		cout << "Solid has " << domain->solids[solid[igroup]]->np
		     << " particles" << endl;
	      }
	    else
	      {
		x    = &domain->solids[solid[igroup]]->grid->x;
		nmax = domain->solids[solid[igroup]]->grid->nnodes_local
		  + domain->solids[solid[igroup]]->grid->nnodes_ghost;
		mask = &domain->solids[solid[igroup]]->grid->mask;
		cout << "Grid has " << domain->solids[solid[igroup]]->grid->nnodes
		     << " nodes" << endl;
	      }

	    int n = 0;

	    for (int ip = 0; ip < nmax; ip++)
	      {
		if (domain->regions[iregion]->match((*x)[ip][0],(*x)[ip][1],(*x)[ip][2])) {
		  (*mask)[ip] |= bit;
		  n++;
		}
	      }

	    int n_tot = 0;
	    MPI_Allreduce(&n,&n_tot,1,MPI_INT,MPI_SUM,universe->uworld);
	
	    cout << n_tot << " " << pon[igroup] << " from solid " << domain->solids[solid[igroup]]->id << " found" << endl;
	  }

      }
    else
      {
	error->all(FLERR, "Error: unknown keyword in group command: " + args[3] + ".\n");
      }
  }
}

/* ----------------------------------------------------------------------
   return group index if name matches existing group, -1 if no such group
------------------------------------------------------------------------- */

int Group::find(string name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (name.compare(names[igroup]) == 0)
      return igroup;
  return -1;
}

/* ----------------------------------------------------------------------
   return index of first available group
   should never be called when group limit has been reached
------------------------------------------------------------------------- */

int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == "")
      return igroup;
  return -1;
}

double Group::xcm(int igroup, int dir)
{
  vector<Eigen::Vector3d> *x;
  vector<double> *mass;
  int nmax;
  vector<int> *mask;
  double com = 0;
  double mass_tot = 0;
  int groupbit    = group->bitmask[igroup];

  if (solid[igroup] == -1)
  {
    // Consider all solids

    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
      {
	if (pon[igroup].compare("particles") == 0)
	  {
	    x    = &domain->solids[solid[igroup]]->x;
	    mass = &domain->solids[solid[igroup]]->mass;
	    nmax = domain->solids[solid[igroup]]->np_local;
	    mask = &domain->solids[solid[igroup]]->mask;
	  }
	else
	  {
	    x    = &domain->solids[solid[igroup]]->grid->x;
	    mass = &domain->solids[solid[igroup]]->grid->mass;
	    nmax = domain->solids[solid[igroup]]->grid->nnodes_local;
	    mask = &domain->solids[solid[igroup]]->grid->mask;
	  }
    
	for (int ip = 0; ip < nmax; ip++)
	  {
	    if ((*mask)[ip] & groupbit)
	      {
		com += (*x)[ip][dir] * (*mass)[ip];
		mass_tot += (*mass)[ip];
	      }
	  }
      }
  }
  else
    {
      int isolid = solid[igroup];

	  if (pon[igroup].compare("particles") == 0)
	    {
	      x    = &domain->solids[solid[igroup]]->x;
	      mass = &domain->solids[solid[igroup]]->mass;
	      nmax = domain->solids[solid[igroup]]->np_local;
	      mask = &domain->solids[solid[igroup]]->mask;
	    }
	  else
	    {
	      x    = &domain->solids[solid[igroup]]->grid->x;
	      mass = &domain->solids[solid[igroup]]->grid->mass;
	      nmax = domain->solids[solid[igroup]]->grid->nnodes_local;
	      mask = &domain->solids[solid[igroup]]->grid->mask;
	    }
    
      for (int ip = 0; ip < nmax; ip++)
	{
	  if ((*mask)[ip] & groupbit)
	    {
	      com += (*x)[ip][dir] * (*mass)[ip];
	      mass_tot += (*mass)[ip];
	    }
	}
    }

  double com_reduced,  mass_tot_reduced;

  MPI_Allreduce(&com,&com_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  MPI_Allreduce(&mass_tot,&mass_tot_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);

  if (mass_tot_reduced) return com_reduced/mass_tot_reduced;
  else return 0;
}

double Group::internal_force(int igroup, int dir)
{
  
  vector<Eigen::Vector3d> *f;
  int nmax;
  vector<int> *mask;
  double resulting_force = 0;
  int groupbit           = group->bitmask[igroup];

  if (solid[igroup] == -1)
  {
    // Consider all solids

    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
      {
	if (pon[igroup].compare("particles") == 0)
	  {
	    f =    &domain->solids[solid[igroup]]->f;
	    nmax = domain->solids[solid[igroup]]->np_local;
	    mask = &domain->solids[solid[igroup]]->mask;
	  }
	else
	  {
	    f =    &domain->solids[solid[igroup]]->grid->f;
	    nmax = domain->solids[solid[igroup]]->grid->nnodes_local;
	    mask = &domain->solids[solid[igroup]]->grid->mask;
	  }
    
	for (int ip = 0; ip < nmax; ip++)
	  {
	    if ((*mask)[ip] & groupbit)
	      {
		resulting_force += (*f)[ip][dir];
	      }
      }
    }
  }
  else
  {
    int isolid = solid[igroup];

    if (pon[igroup].compare("particles") == 0)
      {
	f =    &domain->solids[solid[igroup]]->f;
	nmax = domain->solids[solid[igroup]]->np_local;
	mask = &domain->solids[solid[igroup]]->mask;
      }
    else
      {
	f =    &domain->solids[solid[igroup]]->grid->f;
	nmax = domain->solids[solid[igroup]]->grid->nnodes_local;
	mask = &domain->solids[solid[igroup]]->grid->mask;
      }
    
    for (int ip = 0; ip < nmax; ip++)
      {
	if ((*mask)[ip] & groupbit)
	  {
	    resulting_force += (*f)[ip][dir];
      }
    }
  }

  double resulting_force_reduced;

  MPI_Allreduce(&resulting_force,&resulting_force_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);

  return resulting_force_reduced;
}

double Group::external_force(int igroup, int dir)
{
  if (pon[igroup].compare("nodes") == 0)
    {
      error->all(FLERR, "Error: cannot calculate the external forces applied to the node group "
		 + names[igroup] + ".\n");
    }
  
  vector<Eigen::Vector3d> *f;
  int nmax;
  vector<int> *mask;
  double resulting_force = 0;
  int groupbit           = group->bitmask[igroup];

  if (solid[igroup] == -1)
  {
    // Consider all solids

    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
      {
	f =    &domain->solids[solid[igroup]]->f;
	nmax = domain->solids[solid[igroup]]->np_local;
	mask = &domain->solids[solid[igroup]]->mask;
    
	for (int ip = 0; ip < nmax; ip++)
	  {
	    if ((*mask)[ip] & groupbit)
	      {
		resulting_force += (*f)[ip][dir];
	      }
      }
    }
  }
  else
  {
    int isolid = solid[igroup];

    f =    &domain->solids[solid[igroup]]->f;
    nmax = domain->solids[solid[igroup]]->np_local;
    mask = &domain->solids[solid[igroup]]->mask;
    
    for (int ip = 0; ip < nmax; ip++)
      {
	if ((*mask)[ip] & groupbit)
	  {
	    resulting_force += (*f)[ip][dir];
      }
    }
  }

  double resulting_force_reduced;

  MPI_Allreduce(&resulting_force,&resulting_force_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);

  return resulting_force;
}
