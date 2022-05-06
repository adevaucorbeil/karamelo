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
#include <matrix.h>

class Expression;

class FixInitialVelocityParticles : public Fix {
 public:
  FixInitialVelocityParticles(MPM *, vector<string>);

  void prepare() {};

  void initial_integrate(Solid &solidp);

  void write_restart(ofstream *) {};
  void read_restart(ifstream *) {};

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, initial_velocity_particles, group, vx)\n"},
      {2, "Usage: fix(fix-ID, initial_velocity_particles, group, vx, vy)\n"},
      {3, "Usage: fix(fix-ID, initial_velocity_particles, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};
  
    bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?

  Expression *v[3];
};

#endif
#endif

