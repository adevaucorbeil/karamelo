/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(velocity,FixVelocity)

#else

#ifndef MPM_FIX_VELOCITY_H
#define MPM_FIX_VELOCITY_H

#include "fix.h"
#include <vector>

class FixVelocity : public Fix {
 public:
  FixVelocity(class MPM *, vector<string>);
  ~FixVelocity();
  // int setmask();
  void init();
  void setup();
  // void min_setup(int);
  // void initial_integrate(int);
  // void post_force(int);
  // double compute_vector(int);
  // double memory_usage();

 private:
  // double xvalue,yvalue,zvalue;
  // int varflag,iregion;
  // char *xstr,*ystr,*zstr;
  // char *idregion;
  // int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  // double foriginal[3],foriginal_all[3];
  // int force_flag;
  // int nlevels_respa;

  // int maxatom;
  // double **sforce;
};

#endif
#endif

