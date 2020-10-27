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
#include <Eigen/Eigen>
#include "fix_check_solution.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "output.h"
#include "math_special.h"
#include "universe.h"
#include "solid.h"
#include "error.h"

using namespace std;
using namespace FixConst;
using namespace MathSpecial;
using namespace Eigen;

FixChecksolution::FixChecksolution(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    xset = yset = zset = false;
    return;
  }

  if (domain->dimension == 3 && args.size()<6) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 2 && args.size()<5) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 5 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"_check_solution needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixChecksolution with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  if (args[3].compare("NULL") != 0) {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      yvalue = input->parsev(args[4]);
      yset = true;
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zvalue = input->parsev(args[5]);
      zset = true;
    }
  }
}

FixChecksolution::~FixChecksolution()
{
}

void FixChecksolution::init()
{
}

void FixChecksolution::setup()
{
}

void FixChecksolution::setmask() {
  mask = 0;
  mask |= FINAL_INTEGRATE;
}


void FixChecksolution::final_integrate() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps) return;
  // cout << "In FixChecksolution::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;
    
  int solid = group->solid[igroup];

  Solid *s;

  Eigen::Vector3d error, error_reduced;
  Eigen::Vector3d u_th;

  error.setZero();
  u_th.setZero();

  double vtot;

  if (solid == -1) {
    vtot = 0;
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      vtot += s->vtot;

      for (int in = 0; in < s->np_local; in++) {
	if (s->mask[in] & groupbit) {
	  if (xset) {
	    (*input->vars)["x0"] = Var("x0", s->x0[in][0]);
	    ux = xvalue.result(mpm);
	    error[0] += s->vol0[in]*square(ux-(s->x[in][0]-s->x0[in][0]));
	    u_th[0] += s->vol0[in]*ux*ux;
	  }
	  if (yset) {
	    (*input->vars)["y0"] = Var("y0", s->x0[in][1]);
	    uy = yvalue.result(mpm);
	    error[1] += s->vol0[in]*square(uy-(s->x[in][1]-s->x0[in][1]));
	    u_th[1] += s->vol0[in]*uy*uy;
	  }
	  if (zset) {
	    (*input->vars)["z0"] = Var("z0", s->x0[in][2]);
	    uz = zvalue.result(mpm);
	    error[2] += s->vol0[in]*square(uz-(s->x[in][2]-s->x0[in][2]));
	    u_th[2] += s->vol0[in]*uz*uz;
	  }
	}
      }
    }
  } else {
    s = domain->solids[solid];
    vtot += s->vtot;

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
	if (xset) {
	  (*input->vars)["x0"] = Var("x0", s->x0[in][0]);
	  ux = xvalue.result(mpm);
	  error[0] += s->vol0[in]*square(ux-(s->x[in][0]-s->x0[in][0]));
	  u_th[0] += s->vol0[in]*ux*ux;
	}
	if (yset) {
	  (*input->vars)["y0"] = Var("y0", s->x0[in][1]);
	  uy = yvalue.result(mpm);
	  error[1] += s->vol0[in]*square(uy-(s->x[in][1]-s->x0[in][1]));
	  u_th[1] += s->vol0[in]*uy*uy;
	}
	if (zset) {
	  (*input->vars)["z0"] = Var("z0", s->x0[in][2]);
	  uz = zvalue.result(mpm);
	  error[2] += s->vol0[in]*square(uz-(s->x[in][2]-s->x0[in][2]));
	  u_th[2] += s->vol0[in]*uz*uz;
	}
      }
    }

  }

  // Reduce error:
  MPI_Allreduce(error.data(),error_reduced.data(),3,MPI_DOUBLE,MPI_SUM,universe->uworld);

  (*input->vars)[id+"_s"]=Var(id+"_s", sqrt((error_reduced[0] + error_reduced[1] + error_reduced[2])/vtot));
  (*input->vars)[id+"_x"]=Var(id+"_x", (*input->vars)[id+"_x"].result() + update->dt*(error_reduced[0] + error_reduced[1] + error_reduced[2]));
  (*input->vars)[id+"_y"]=Var(id+"_y", (*input->vars)[id+"_y"].result() + update->dt*(u_th[0] + u_th[1] + u_th[2]));
  (*input->vars)[id+"_z"]=Var(id+"_z", sqrt((*input->vars)[id+"_x"].result()/(*input->vars)[id+"_y"].result()));
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}


void FixChecksolution::write_restart(ofstream *of) {
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

void FixChecksolution::read_restart(ifstream *ifr) {
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
