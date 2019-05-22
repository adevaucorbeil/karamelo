/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef SOLID_CLASS

SolidStyle(rectangle,SolRectangle)

#else

#ifndef MPM_SOLID_RECTANGLE_H
#define MPM_SOLID_RECTANGLE_H

#include "solid.h"

class SolRectangle : public Solid {

 public:
  SolRectangle(class MPM *, vector<string>);
  ~SolRectangle();

};

#endif
#endif
