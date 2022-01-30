/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(temperature_particles,FixTemperatureParticles)

#else

#ifndef MPM_FIX_TEMPERATURE_PARTICLES_H
#define MPM_FIX_TEMPERATURE_PARTICLES_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixTemperatureParticles : public Fix {
 public:
  FixTemperatureParticles(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void initial_integrate(Solid &solid, int ip);
  void post_advance_particles(Solid &solid, int ip);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:  
  string usage = "Usage: fix(fix-ID, temperature_nodes, group-ID, T)\n";
  int Nargs = 4;

  Var Tvalue;                      //< Temperature variable.
  Var Tprevvalue;                  //< Temperature variable from previous time step.
};

#endif
#endif

