/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_EOS_H
#define MPM_EOS_H

#include "pointers.h"
#include <vector>
#include <Eigen/Eigen>

class EOS : protected Pointers {
 public:
  string id;

  EOS(class MPM *, vector<string>);
  virtual ~EOS();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each EOS
  //virtual compute_pressure()
  virtual double rho0() = 0;
  virtual double K() = 0;
  virtual double compute_pressure(const double J, const double rho, const double e, const double damage) = 0;
  //protected:
};

#endif
