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

#include <matrix.h>
#include <string>
#include <group.h>
#include <domain.h>
#include <region.h>
#include <solid.h>
#include <error.h>
#include <universe.h>

#define MAX_GROUP 32

using namespace std;


Group::Group(MPM *mpm) : Pointers(mpm)
{
  names       = new string[MAX_GROUP];
  bitmask     = new int[MAX_GROUP];
  //inversemask = new int[MAX_GROUP];
  pon         = new string[MAX_GROUP];
  solid       = new int[MAX_GROUP];
  region      = new int[MAX_GROUP];
  n_tot_      = new int[MAX_GROUP];

  for (int i = 0; i < MAX_GROUP; i++) {
    names[i] = "";
    bitmask[i] = 1 << i;
    //inversemask[i] = bitmask[i] ^ ~0;
    pon[i] = "all";
    solid[i] = -1;
    region[i] = -1;
    n_tot_[i] = 0;
  }

  // create "all" group
  names[0] = "all";
  ngroup   = 1;
}

Group::~Group()
{
  delete[] names;
  delete[] bitmask;
  //delete[] inversemask;
  delete[] pon;
  delete[] solid;
  delete[] region;
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
    if (ngroup == MAX_GROUP) {
      error->all(FLERR, "Too many groups.\n");
    }
    igroup        = find_unused();
    names[igroup] = args[0];
    ngroup++;
  }

  int bit = bitmask[igroup];

  // Operates on particles or nodes:
  if (args[1] == "particles")
  {
    pon[igroup] = "particles";
  }
  else if (args[1] == "nodes")
  {
    pon[igroup] = "nodes";
  } else {
    error->all(FLERR,"Error: do not understand keyword " + args[1] + ", \"particles\" or \"nodes\" expected.\n");
  }

  // style = region
  // add to group if atom is in region

  if (args[2] == "region")
  {
    // Look for the region ID (if exists):
    region[igroup] = domain->find_region(args[3]);
    if (region[igroup] == -1)
      {
      error->all(FLERR, "Error: could not find region " + args[3] + ".\n");
      }

    if (args[4] == "all")
      {

	/* For all particles of all solids, check if they are in the region.
	   If so asign them the right mask */
	solid[igroup] = -1; // Consider all solids

	for (int isolid = 0; isolid < domain->solids.size(); isolid++)
	  {

	    vector<Vector3d> *x;
	    int nmax;
	    vector<int> *mask;

	    if (pon[igroup] == "particles")
	      {
		x    = &domain->solids[isolid]->x0;
		nmax = domain->solids[isolid]->np_local;
		mask = &domain->solids[isolid]->mask;

                // if (universe->me == 0) {
                //   cout << "Solid has " << domain->solids[isolid]->np
                //        << " particles" << endl;
                // }
              }
	    else
	      {
		x    = &domain->solids[isolid]->grid->x0;
		nmax = domain->solids[isolid]->grid->nnodes_local
		  + domain->solids[isolid]->grid->nnodes_ghost;
		mask = &domain->solids[isolid]->grid->mask;

                // if (universe->me == 0) {
                //   cout << "Grid has " << domain->solids[isolid]->grid->nnodes
                //        << " nodes" << endl;
                // }
              }

	    int n = 0;

	    for (int ip = 0; ip < nmax; ip++)
	      {
		if (domain->regions[region[igroup]]->match((*x)[ip][0],(*x)[ip][1],(*x)[ip][2]))
		  {
		    (*mask)[ip] |= bit;
		    n++;
		  }
	      }

	    n_tot_[igroup] = 0;
	    MPI_Allreduce(&n,&n_tot_[igroup],1,MPI_INT,MPI_SUM,universe->uworld);

            if (universe->me == 0) {
              cout << n_tot_[igroup] << " " << pon[igroup] << " from solid "
                   << domain->solids[isolid]->id << " found" << endl;
            }
          }
      }
    else if (args[4] == "solid")
      {

	for (int i=5; i<args.size(); i++)
	  {
	    solid[igroup] = domain->find_solid(args[i]);
	    if (solid[igroup] == -1) {
	      error->all(FLERR, "Error: cannot find solid with ID " + args[i] + ".\n");
	    }

	    vector<Vector3d> *x;
	    int nmax;
	    vector<int> *mask;

	    if (pon[igroup] == "particles")
	      {
		x    = &domain->solids[solid[igroup]]->x0;
		nmax = domain->solids[solid[igroup]]->np_local;
		mask = &domain->solids[solid[igroup]]->mask;
		// if (universe->me == 0) {
		//   cout << "Solid has " << domain->solids[solid[igroup]]->np
		//        << " particles" << endl;
		// }
	      }
	    else
	      {
		x    = &domain->solids[solid[igroup]]->grid->x0;
		nmax = domain->solids[solid[igroup]]->grid->nnodes_local
		  + domain->solids[solid[igroup]]->grid->nnodes_ghost;
		mask = &domain->solids[solid[igroup]]->grid->mask;
		// if (universe->me == 0) {
		//   cout << "Grid has " << domain->solids[solid[igroup]]->grid->nnodes
		//        << " nodes" << endl;
		// }
	      }

	    int n = 0;

	    for (int ip = 0; ip < nmax; ip++)
	      {
		if (domain->regions[region[igroup]]->match((*x)[ip][0],(*x)[ip][1],(*x)[ip][2])) {
		  (*mask)[ip] |= bit;
		  n++;
		}
	      }

	    n_tot_[igroup] = 0;
	    MPI_Allreduce(&n,&n_tot_[igroup],1,MPI_INT,MPI_SUM,universe->uworld);

            if (universe->me == 0) {
              cout << n_tot_[igroup] << " " << pon[igroup] << " from solid "
                   << domain->solids[solid[igroup]]->id << " found" << endl;
            }
          }

      }
    else
      {
	error->all(FLERR, "Error: unknown keyword in group command: " + args[3] + ".\n");
      }
  } else {
    error->all(FLERR,
               "Error: unknown keyword in group command: " + args[2] + ".\n");
  }
}

/* ----------------------------------------------------------------------
   return group index if name matches existing group, -1 if no such group
------------------------------------------------------------------------- */

int Group::find(string name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (name == names[igroup])
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
  vector<Vector3d> *x;
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
	if (pon[igroup] == "particles")
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

	  if (pon[igroup] == "particles")
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
  
  vector<Vector3d> *f;
  int nmax;
  vector<int> *mask;
  double resulting_force = 0;
  int groupbit           = group->bitmask[igroup];

  if (solid[igroup] == -1)
  {
    // Consider all solids

    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
      {
	if (pon[igroup] == "particles")
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

    if (pon[igroup] == "particles")
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
  if (pon[igroup] == "nodes")
    {
      error->all(FLERR, "Error: cannot calculate the external forces applied to the node group "
		 + names[igroup] + ".\n");
    }
  
  vector<Vector3d> *f;
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


/*! Write groups to restart file
 */
void Group::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&ngroup), sizeof(int));

  if (ngroup <= 1) return;

  size_t Nr = 0;
  bool p_or_n = false;

  for (int igroup = 1; igroup < ngroup; igroup++) {
    Nr = names[igroup].size();
    of->write(reinterpret_cast<const char *>(&Nr), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(names[igroup].c_str()), Nr);
    // cout << "Group name = " << names[igroup] << endl;

    of->write(reinterpret_cast<const char *>(&bitmask[igroup]), sizeof(int));
    //of->write(reinterpret_cast<const char *>(&inversemask[igroup]), sizeof(int));

    if (pon[igroup] == "particles")
      p_or_n = false;
    else 
      p_or_n = true;

    if (pon[igroup] == "all") {
      error->all(FLERR,"Error: pon==all \n");
    }

    of->write(reinterpret_cast<const char *>(&p_or_n), sizeof(bool));

    of->write(reinterpret_cast<const char *>(&solid[igroup]), sizeof(int));
    // cout << "Group solid = " << solid[igroup] << endl;
    of->write(reinterpret_cast<const char *>(&region[igroup]), sizeof(int));
    // cout << "Group region = " << region[igroup] << endl;
    // cout << "Group pon = " << pon[igroup] << endl;
  }
}

/*! Read groups from restart file
 */
void Group::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&ngroup), sizeof(int));

  if (ngroup <= 1) return;

  size_t Nr = 0;
  bool p_or_n = false;

  for (int igroup = 1; igroup < ngroup; igroup++) {
    ifr->read(reinterpret_cast<char *>(&Nr), sizeof(size_t));
    names[igroup].resize(Nr);

    ifr->read(reinterpret_cast<char *>(&names[igroup][0]), Nr);
    // cout << "Group name = " << names[igroup] << endl;

    ifr->read(reinterpret_cast<char *>(&bitmask[igroup]), sizeof(int));
    ifr->read(reinterpret_cast<char *>(&p_or_n), sizeof(bool));

    if (p_or_n == false) {
      pon[igroup] = "particles";
    } else {
      pon[igroup] = "nodes";
    }

    ifr->read(reinterpret_cast<char *>(&solid[igroup]), sizeof(int));
    ifr->read(reinterpret_cast<char *>(&region[igroup]), sizeof(int));

    // cout << "Group solid = " << solid[igroup] << endl;
    // cout << "Group region = " << region[igroup] << endl;
    // cout << "Group pon = " << pon[igroup] << endl;
    if (region[igroup] == -1) {
      error->all(FLERR, "Error: could not find region with ID " + to_string(region[igroup]) + ".\n");
    }

    if (solid[igroup] == -1) {
      // Consider all solids
      for (int isolid = 0; isolid < domain->solids.size(); isolid++) {

        vector<Vector3d> *x;
        int nmax;
        vector<int> *mask;

        if (pon[igroup] == "particles") {
          x = &domain->solids[isolid]->x0;
          nmax = domain->solids[isolid]->np_local;
          mask = &domain->solids[isolid]->mask;
	  if (universe->me == 0) {
	    cout << "Solid has " << domain->solids[isolid]->np << " particles"
		 << endl;
	  }
        } else {
          x = &domain->solids[isolid]->grid->x0;
          nmax = domain->solids[isolid]->grid->nnodes_local +
                 domain->solids[isolid]->grid->nnodes_ghost;
          mask = &domain->solids[isolid]->grid->mask;
	  if (universe->me == 0) {
	    cout << "Grid has " << domain->solids[isolid]->grid->nnodes
		 << " nodes" << endl;
	  }
        }

        int n = 0;

        for (int ip = 0; ip < nmax; ip++) {
          if (domain->regions[region[igroup]]->match((*x)[ip][0], (*x)[ip][1],
                                              (*x)[ip][2])) {
            (*mask)[ip] |= bitmask[igroup];
            n++;
          }
        }

        n_tot_[igroup] = 0;
        MPI_Allreduce(&n, &n_tot_[igroup], 1, MPI_INT, MPI_SUM, universe->uworld);

	if (universe->me == 0) {
	  cout << n_tot_[igroup] << " " << pon[igroup] << " from solid "
	       << domain->solids[isolid]->id << " found" << endl;
	}
      }
    } else {
      vector<Vector3d> *x;
      int nmax;
      vector<int> *mask;

      if (pon[igroup] == "particles") {
        x = &domain->solids[solid[igroup]]->x0;
        nmax = domain->solids[solid[igroup]]->np_local;
        mask = &domain->solids[solid[igroup]]->mask;
	if (universe->me == 0) {
	  cout << "Solid has " << domain->solids[solid[igroup]]->np
	       << " particles" << endl;
	}
      } else {
        x = &domain->solids[solid[igroup]]->grid->x0;
        nmax = domain->solids[solid[igroup]]->grid->nnodes_local +
               domain->solids[solid[igroup]]->grid->nnodes_ghost;
        mask = &domain->solids[solid[igroup]]->grid->mask;
	if (universe->me == 0) {
	  cout << "Grid has " << domain->solids[solid[igroup]]->grid->nnodes
	       << " nodes" << endl;
	}
      }

      int n = 0;

      for (int ip = 0; ip < nmax; ip++) {
        if (domain->regions[region[igroup]]->match((*x)[ip][0], (*x)[ip][1],
                                            (*x)[ip][2])) {
          (*mask)[ip] |= bitmask[igroup];
          n++;
        }
      }

      n_tot_[igroup] = 0;
      MPI_Allreduce(&n, &n_tot_[igroup], 1, MPI_INT, MPI_SUM, universe->uworld);

      if (universe->me == 0) {
	cout << n_tot_[igroup] << " " << pon[igroup] << " from solid "
	     << domain->solids[solid[igroup]]->id << " found" << endl;
      }
    }
  }
}
