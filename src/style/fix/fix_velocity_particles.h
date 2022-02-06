///* -*- c++ -*- ----------------------------------------------------------*/
//
//#ifdef FIX_CLASS
//
//FixStyle(velocity_particles,FixVelocityParticles)
//
//#else
//
//#ifndef MPM_FIX_VELOCITY_PARTICLES_H
//#define MPM_FIX_VELOCITY_PARTICLES_H
//
//#include <fix.h>
//#include <var.h>
//#include <matrix.h>
//
//class FixVelocityParticles : public Fix {
// public:
//  FixVelocityParticles(MPM *, vector<string>);
//
//  void prepare();
//  void reduce();
//  
//  void initial_integrate(Solid &solid, int ip);
//  void post_advance_particles(Solid &solid, int ip);
//
//  void write_restart(ofstream *);
//  void read_restart(ifstream *);
//
//private:
//  const map<int, string> usage = {
//      {1, "Usage: fix(fix-ID, velocity_particles, group, vx)\n"},
//      {2, "Usage: fix(fix-ID, velocity_particles, group, vx, vy)\n"},
//      {3, "Usage: fix(fix-ID, velocity_particles, group, vx, vy, vz)\n"}};
//  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};
//
//  Var xvalue, yvalue, zvalue;                  //< Velocities in x, y, and z directions.
//  Var xprevvalue, yprevvalue, zprevvalue;      //< Velocities in x, y, and z directions from previous time step.
//  bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?
//
//  vector<Vector3d> xold;          // particles' old position
//  Vector3d ftot;
//};
//
//#endif
//#endif
//
