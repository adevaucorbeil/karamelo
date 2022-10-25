/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(temperature_particles,FixTemperatureParticles)

#else

#ifndef MPM_FIX_TEMPERATURE_PARTICLES_H
#define MPM_FIX_TEMPERATURE_PARTICLES_H

#include <fix.h>
#include <matrix.h>

class Expression;

class FixTemperatureParticles : public Fix {
 public:
  FixTemperatureParticles(MPM *, vector<string>);

  void prepare() {};
  void reduce() {};
  
  void initial_integrate(Solid &solid);
  void post_advance_particles(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:  
  string usage = "Usage: fix(fix-ID, temperature_nodes, group-ID, T)\n";
  int Nargs = 4;

  Expression *Tvalue;                      //< Temperature variable.
  Expression *Tprevvalue;                  //< Temperature variable from previous time step.
};

#endif
#endif
