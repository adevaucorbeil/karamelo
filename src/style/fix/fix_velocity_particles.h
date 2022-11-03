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

FixStyle(velocity_particles,FixVelocityParticles)

#else

#ifndef MPM_FIX_VELOCITY_PARTICLES_H
#define MPM_FIX_VELOCITY_PARTICLES_H

#include <fix.h>
#include <matrix.h>

class Expression;

class FixVelocityParticles : public Fix {
public:
  FixVelocityParticles(MPM *, vector<string>);

  void prepare();
  void reduce();
 
  void initial_integrate(Solid &solid);
  void post_advance_particles(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
    {1, "Usage: fix(fix-ID, velocity_particles, group, vx)\n"},
    {2, "Usage: fix(fix-ID, velocity_particles, group, vx, vy)\n"},
    {3, "Usage: fix(fix-ID, velocity_particles, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  vector<Vector3d> xold;          // particles' old position
  Expression *v[3];
  Expression *v_prev[3];
  Vector3d ftot;
};

#endif
#endif

