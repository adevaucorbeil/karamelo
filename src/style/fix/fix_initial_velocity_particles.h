/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef FIX_CLASS

FixStyle(initial_velocity_particles,FixInitialVelocityParticles)

#else

#ifndef MPM_FIX_INITIAL_VELOCITY_PARTICLES_H
#define MPM_FIX_INITIAL_VELOCITY_PARTICLES_H

#include <fix.h>
#include <var.h>
#include <vector>

class FixInitialVelocityParticles : public Fix {
 public:
  FixInitialVelocityParticles(MPM *, vector<string>);
  
  void initial_integrate();

  void write_restart(ofstream *) {};
  void read_restart(ifstream *) {};

private:
  string usage = "Usage: fix(fix-ID, initial_velocity_particles, group-ID, vx, vy, vz)\n";
  int Nargs = 6;
  Var xvalue, yvalue, zvalue;    // Set velocities in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z velocities of the group?
};

#endif
#endif

