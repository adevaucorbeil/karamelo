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

#include <iostream>
#include <vector>
#include <string>
#include <matrix.h>
#include <fix_force_nodes.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <input.h>
#include <universe.h>
#include <grid.h>
#include <error.h>

using namespace std;
using namespace FixConst;


FixForceNodes::FixForceNodes(MPM *mpm, vector<string> args):
  Fix(mpm, args, POST_PARTICLES_TO_GRID)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];

    xset = yset = zset = false;
    return;
  }

  if (args.size() < 6) {
    error->all(FLERR,"Error: too few arguments for fix_force_nodes: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->all(FLERR,"fix_force_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixForceNodes with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (args[4] != "NULL") {
    yvalue = input->parsev(args[4]);
    yset = true;
  }

  if (args[5] != "NULL") {
    zvalue = input->parsev(args[5]);
    zset = true;
  }
}

void FixForceNodes::post_particles_to_grid() {
  // cout << "In FixForceNodes::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double fx, fy, fz;

  if (xset) {
    fx = xvalue.result(mpm);
    // cout << "Set v_update[0] to " << xvalue.eq() << "=" << fx << endl;
  }

  if (yset) {
    fy = yvalue.result(mpm);
    // cout << "Set v_update[1] to " << "=" <<  fy << endl;
  }

  if (zset) {
    fz = zvalue.result(mpm);
    // cout << "Set v_update[2] to " << "=" <<  fz << endl;
  }
    
  int solid = group->solid[igroup];
  Grid *g;

  int n = 0;
  Vector3d ftot, ftot_reduced;
  ftot = Vector3d();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;
      n = 0;

      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
	if (g->mass[in] > 0) {
	  if (g->mask[in] & groupbit) {
	    n++;
	  }
	}
      }
	    
      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
	if (g->mass[in] > 0) {
	  if (g->mask[in] & groupbit) {
	    if (xset) {
	      g->mb[in][0] += fx/((double) n);
	      if(in < g->nnodes_local) ftot[0] += fx/((double) n);
	    }
	    if (yset) {
	      g->mb[in][1] += fy/((double) n);
	      if(in < g->nnodes_local) ftot[1] += fy/((double) n);
	    }
	    if (zset) {
	      g->mb[in][2] += fz/((double) n);
	      if(in < g->nnodes_local) ftot[2] += fz/((double) n);
	    }
	  }
	}
      }
    }
  } else {

    g = domain->solids[solid]->grid;
    n = 0;

    for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
      if (g->mass[in] > 0) {
	if (g->mask[in] & groupbit) {
	  n++;
	}
      }
    }

    for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
      if (g->mass[in] > 0) {
	if (g->mask[in] & groupbit) {
	  if (xset) {
	    g->mb[in][0] += fx/((double) n);
	    if(in < g->nnodes_local) ftot[0] += fx/((double) n);
	  }
	  if (yset) {
	    g->mb[in][1] += fy/((double) n);
	    if(in < g->nnodes_local) ftot[1] += fy/((double) n);
	  }
	  if (zset) {
	    g->mb[in][2] += fz/((double) n);
	    if(in < g->nnodes_local) ftot[2] += fz/((double) n);
	  }
	}
      }
    }
  }

  // Reduce ftot:
  MPI_Allreduce(ftot.elements,ftot_reduced.elements,3,MPI_DOUBLE,MPI_SUM,universe->uworld);

  if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot_reduced[0]);
  if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot_reduced[1]);
  if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot_reduced[2]);
  // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}

void FixForceNodes::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset) {
    string eq = xvalue.eq();
    size_t N = eq.size();
    double value = xvalue.result();
    bool cst = xvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }

  if (yset) {
    string eq = yvalue.eq();
    size_t N = eq.size();
    double value = yvalue.result();
    bool cst = yvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }

  if (zset) {
    string eq = zvalue.eq();
    size_t N = eq.size();
    double value = zvalue.result();
    bool cst = zvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }
}

void FixForceNodes::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    xvalue = Var(eq, value, cst);
  }

  if (yset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    yvalue = Var(eq, value, cst);
  }

  if (zset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    zvalue = Var(eq, value, cst);
  }
}
