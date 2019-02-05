/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include "pointers.h"
#include "eos.h"
#include "grid.h"
#include <vector>

class Solid : protected Pointers {
 public:
  string id;                          // solid id
  vector< array<double, 3> > x;       // particles' current position
  vector< array<double, 3> > x0;      // particles' reference position
  vector< double > vol0;              // particles' reference volume
  vector< double > vol;               // particles' current volume
  vector< double > mass;              // particles' current mass

  class EOS *eos;                     // Equation-of-State

  class Grid *grid;                    // background grid

  Solid(class MPM *, vector<string>);
  virtual ~Solid();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);
  void grow(int);

protected:
  bigint nparticles; // number of particles
};

#endif
