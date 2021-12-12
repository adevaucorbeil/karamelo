/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2020) Alban de Vaucorbeil, alban.devaucorbeil@deakin.edu.au
 * Institute for Frontier Materials, Deakin University
 * Geelong VIC 3216, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <fix_velocity_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;


FixVelocityParticles::FixVelocityParticles(MPM *mpm, vector<string> args) : Fix(mpm, args)
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

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_velocity_nodes.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_velocity_nodes needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixVelocityParticles with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;


  string time = "time";


  if (args[3].compare("NULL") != 0) {
    // xvalue = input->parsev(args[3]);
    xset = true;
    xvalue = input->parsev(args[3]);

    string previous = args[3];

    // Replace "time" by "time - dt" in the x argument:
    previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
    xprevvalue = input->parsev(previous);
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      yvalue = input->parsev(args[4]);
      yset = true;

      string previous = args[4];

      // Replace "time" by "time - dt" in the y argument:
      previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
      yprevvalue = input->parsev(previous);
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zvalue = input->parsev(args[5]);
      zset = true;

      string previous = args[5];

      // Replace "time" by "time - dt" in the z argument:
      previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
      zprevvalue = input->parsev(previous);
    }
  }
}

FixVelocityParticles::~FixVelocityParticles()
{
}

void FixVelocityParticles::init()
{
}

void FixVelocityParticles::setup()
{
}

void FixVelocityParticles::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_ADVANCE_PARTICLES;
}


void FixVelocityParticles::initial_integrate() {
  // Go through all the particles in the group and set v_update to the right value:
  Vector3d xtemp;

  int solid = group->solid[igroup];
  Solid *s;

  int n = 0;
  xold.clear();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      n = 0;

      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mask[ip] & groupbit) {
	  xtemp = s->x[ip];
	  (*input->vars)["x"] = Var("x", s->x[ip][0]);
	  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	  (*input->vars)["y"] = Var("y", s->x[ip][1]);
	  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	  (*input->vars)["z"] = Var("z", s->x[ip][2]);
	  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

	  if (xset) {
	    s->v_update[ip][0] = xvalue.result(mpm);
	    s->v[ip][0] = xprevvalue.result(mpm);
	  }
	  if (yset) {
	    s->v_update[ip][1] = yvalue.result(mpm);
	    s->v[ip][1] = yprevvalue.result(mpm);
	  }
	  if (zset) {
	    s->v_update[ip][2] = zvalue.result(mpm);
	    s->v[ip][2] = zprevvalue.result(mpm);
	  }
	  // if (s->ptag[ip] == 4371) {
	  //   printf("fix: v=[%4.3e %4.3e %4.3e]\tv_update=[%4.3e %4.3e %4.3e]\ta=[%4.3e %4.3e %4.3e]\n", s->v[ip][0], s->v[ip][1], s->v[ip][2], s->v_update[ip][0], s->v_update[ip][1], s->v_update[ip][2], s->a[ip][0], s->a[ip][1], s->a[ip][2]);
	  // }
	  xold.push_back(xtemp);
	  n++;
	}
      }
      // cout << "v_update for " << n << " particles from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
	xtemp = s->x[ip];
	(*input->vars)["x"] = Var("x", s->x[ip][0]);
	(*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	(*input->vars)["y"] = Var("y", s->x[ip][1]);
	(*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	(*input->vars)["z"] = Var("z", s->x[ip][2]);
	(*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

	if (xset) {
	  s->v_update[ip][0] = xvalue.result(mpm);
	  s->v[ip][0] = xprevvalue.result(mpm);
	}
	if (yset) {
	  s->v_update[ip][1] = yvalue.result(mpm);
	  s->v[ip][1] = yprevvalue.result(mpm);
	}
	if (zset) {
	  s->v_update[ip][2] = zvalue.result(mpm);
	  s->v[ip][2] = zprevvalue.result(mpm);
	}
	// if (s->ptag[ip] == 4371) {
	//   printf("fix: v=[%4.3e %4.3e %4.3e]\tv_update=[%4.3e %4.3e %4.3e]\ta=[%4.3e %4.3e %4.3e]\n", s->v[ip][0], s->v[ip][1], s->v[ip][2], s->v_update[ip][0], s->v_update[ip][1], s->v_update[ip][2], s->a[ip][0], s->a[ip][1], s->a[ip][2]);
	// }
	xold.push_back(xtemp);
	n++;
      }
    }
    // cout << "v_update for " << n << " particles from solid " << domain->solids[solid]->id << " set." << endl;
  }
}

void FixVelocityParticles::post_advance_particles() {
  // Go through all the particles in the group and set v to the right value:
  double vx, vy, vz;
  
  int solid = group->solid[igroup];
  Solid *s;
  Eigen::Vector3d Dv, ftot, ftot_reduced;

  int n = 0;
  ftot.setZero();
  double inv_dt = 1.0/update->dt;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      n = 0;

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mask[ip] & groupbit) {
          Dv.setZero();
	  (*input->vars)["x"] = Var("x", xold[n][0]);
	  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	  (*input->vars)["y"] = Var("y", xold[n][1]);
	  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	  (*input->vars)["z"] = Var("z", xold[n][2]);
	  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

          if (xset) {
            vx = xvalue.result(mpm);
            Dv[0] = vx - s->v[ip][0];
            s->v[ip][0] = vx;
            s->x[ip][0] = xold[n][0] + update->dt * vx;
          }
          if (yset) {
            vy = yvalue.result(mpm);
            Dv[1] = vy - s->v[ip][1];
            s->v[ip][1] = vy;
            s->x[ip][1] = xold[n][1] + update->dt * vy;
          }
          if (zset) {
            vz = zvalue.result(mpm);
            Dv[2] = vz - s->v[ip][2];
            s->v[ip][2] = vz;
            s->x[ip][2] = xold[n][2] + update->dt * vz;
          }
          ftot += (inv_dt * s->mass[ip]) * Dv;
          n++;
        }
      }
      // cout << "v for " << n << " particles from solid " <<
      // domain->solids[isolid]->id << " set." << endl;
    }
  } else {
    s = domain->solids[solid];
    n = 0;
    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
        Dv.setZero();
	(*input->vars)["x"] = Var("x", xold[n][0]);
	(*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	(*input->vars)["y"] = Var("y", xold[n][1]);
	(*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	(*input->vars)["z"] = Var("z", xold[n][2]);
	(*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

        if (xset) {
          vx = xvalue.result(mpm);
          Dv[0] = vx - s->v[ip][0];
          s->v[ip][0] = vx;
          s->x[ip][0] = xold[n][0] + update->dt * vx;
        }
        if (yset) {
          vy = yvalue.result(mpm);
          Dv[1] = vy - s->v[ip][1];
          s->v[ip][1] = vy;
          s->x[ip][1] = xold[n][1] + update->dt * vy;
        }
        if (zset) {
          vz = zvalue.result(mpm);
          Dv[2] = vz - s->v[ip][2];
          s->v[ip][2] = vz;
          s->x[ip][2] = xold[n][2] + update->dt * vz;
        }
        ftot += (inv_dt * s->mass[ip]) * Dv;
        n++;
      }
    }
    // cout << "v for " << n << " particles from solid " <<
    // domain->solids[solid]->id << " set." << endl;
  }

  // Reduce ftot:
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixVelocityParticles::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.write_to_restart(of);
    xprevvalue.write_to_restart(of);
  }
  if (yset) {
    yvalue.write_to_restart(of);
    yprevvalue.write_to_restart(of);
  }
  if (zset) {
    zvalue.write_to_restart(of);
    zprevvalue.write_to_restart(of);
  }
}

void FixVelocityParticles::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.read_from_restart(ifr);
    xprevvalue.read_from_restart(ifr);
  }
  if (yset) {
    yvalue.read_from_restart(ifr);
    yprevvalue.read_from_restart(ifr);
  }
  if (zset) {
    zvalue.read_from_restart(ifr);
    zprevvalue.read_from_restart(ifr);
  }
}
