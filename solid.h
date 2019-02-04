/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include "pointers.h"
#include <vector>

class Solid : protected Pointers {
 public:
  string id;
  vector< array<double, 3> > x;
  vector< array<double, 3> > x0;
  vector< double > vol0;
  vector< double > vol;
  vector< double > mass;

  Solid(class MPM *, vector<string>);
  virtual ~Solid();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);
  void grow(int);

protected:
  bigint nparticles; // number of particles
};

#endif
