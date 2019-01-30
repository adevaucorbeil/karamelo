#include "mpm.h"
#include "solid.h"

Solid::Solid(MPM *mpm) : Pointers(mpm)
{
  np = 0; // number of particles
}


Solid::~Solid()
{
}
