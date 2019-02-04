/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef SOLID_CLASS

SolidStyle(block,SolBlock)

#else

#ifndef MPM_SOLID_BLOCK_H
#define MPM_SOLID_BLOCK_H

#include "solid.h"

class SolBlock : public Solid {

 public:
  SolBlock(class MPM *, vector<string>);
  ~SolBlock();

};

#endif
#endif
