#include "fix_velocity_particles.h"
#include "domain.h"
#include "group.h"
#include "input.h"
#include "special_functions.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;


FixVelocityParticles::FixVelocityParticles(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (domain->dimension == 3 && args.size()<6) {
    cout << "Error: too few arguments for fix_velocity_particles: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 2 && args.size()<5) {
    cout << "Error: too few arguments for fix_velocity_particles: requires at least 5 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 1 && args.size()<4) {
    cout << "Error: too few arguments for fix_velocity_particles: requires at least 4 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    cout << "fix_velocity_particles needs to be given a group of particles" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixVelocityParticles with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  args_previous_step = args;

  string time = "time";

  if (args[3].compare("NULL") != 0) {
    // xvalue = input->parsev(args[3]);
    xpos = 3;
    xset = true;

    // Replace "time" by "time - dt" in the x argument:
    args_previous_step[xpos] = SpecialFunc::replace_all(input->parsev(args_previous_step[xpos]).str(), "time", "(time - dt)");
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      ypos = 4;
      // yvalue = input->parsev(args[4]);
      yset = true;
      // Replace "time" by "time - dt" in the y argument:
      args_previous_step[ypos] = SpecialFunc::replace_all(input->parsev(args_previous_step[ypos]).str(), "time", "(time - dt)");
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zpos = 5;
      // zvalue = input->parsev(args[5]);
      zset = true;

      // Replace "time" by "time - dt" in the z argument:
      args_previous_step[zpos] = SpecialFunc::replace_all(input->parsev(args_previous_step[zpos]).str(), "time", "(time - dt)");
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
  double vx, vy, vz;
  double vx_old, vy_old, vz_old;
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
	    vx = input->parsev(args[xpos]).result(mpm);
	    vx_old = input->parsev(args_previous_step[xpos]).result(mpm);
	    s->v_update[ip][0] = vx;
	    s->v[ip][0] = vx_old;
	  }
	  if (yset) {
	    vy = input->parsev(args[ypos]).result(mpm);
	    vy_old = input->parsev(args_previous_step[ypos]).result(mpm);
	    s->v_update[ip][1] = vy;
	    s->v[ip][1] = vy_old;
	  }
	  if (zset) {
	    vz = input->parsev(args[zpos]).result(mpm);
	    vz_old = input->parsev(args_previous_step[zpos]).result(mpm);
	    s->v_update[ip][2] = vz;
	    s->v[ip][2] = vz_old;
	  }
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
	  vx = input->parsev(args[xpos]).result(mpm);
	  vx_old = input->parsev(args_previous_step[xpos]).result(mpm);
	  s->v_update[ip][0] = vx;
	  s->v[ip][0] = vx_old;
	}
	if (yset) {
	  vy = input->parsev(args[ypos]).result(mpm);
	  vy_old = input->parsev(args_previous_step[ypos]).result(mpm);
	  s->v_update[ip][1] = vy;
	  s->v[ip][1] = vy_old;
	}
	if (zset) {
	  vz = input->parsev(args[zpos]).result(mpm);
	  vz_old = input->parsev(args_previous_step[zpos]).result(mpm);
	  s->v_update[ip][2] = vz;
	  s->v[ip][2] = vz_old;
	}
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
            vx = input->parsev(args[xpos]).result(mpm);
            Dv[0] = vx - s->v[ip][0];
            s->v[ip][0] = vx;
            s->x[ip][0] = xold[n][0] + update->dt * vx;
          }
          if (yset) {
            vy = input->parsev(args[ypos]).result(mpm);
            Dv[1] = vy - s->v[ip][1];
            s->v[ip][1] = vy;
            s->x[ip][1] = xold[n][1] + update->dt * vy;
          }
          if (zset) {
            vz = input->parsev(args[zpos]).result(mpm);
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
          vx = input->parsev(args[xpos]).result(mpm);
          Dv[0] = vx - s->v[ip][0];
          s->v[ip][0] = vx;
          s->x[ip][0] = xold[n][0] + update->dt * vx;
        }
        if (yset) {
          vy = input->parsev(args[ypos]).result(mpm);
          Dv[1] = vy - s->v[ip][1];
          s->v[ip][1] = vy;
          s->x[ip][1] = xold[n][1] + update->dt * vy;
        }
        if (zset) {
          vz = input->parsev(args[zpos]).result(mpm);
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
