/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef SOLID_CLASS

SolidStyle(cylinder,SolCylinder)

#else

#ifndef MPM_SOLID_CYLINDER_H
#define MPM_SOLID_CYLINDER_H

#include "solid.h"

class SolCylinder : public Solid {

 public:
  SolCylinder(class MPM *, vector<string>);
  ~SolCylinder();

};

#endif
#endif
