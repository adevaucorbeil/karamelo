/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_VAR_H
#define MPM_VAR_H

#include "input.h"
#include <vector>

class Var : protected Input {
 public:
  // functions
  Var();
  Var(class MPM *);
  Var(class MPM*, double, bool c = true);
  Var(class MPM*, string, double, bool c = false);
  ~Var() {};

  void evaluate();
  double result();
  double result() const {return value;};
  string str() const;
  string eq() const {return equation;};
  bool is_constant() const {return constant;};
  MPM* mpm() const;
  Var operator+(const Var&);
  Var operator-(const Var&);
  Var operator-();
  Var operator*(const Var&);
  Var operator/(const Var&);
  Var operator^(const Var&);
  Var operator=(const Var&);

protected:
  class MPM* mpm_ptr;          // Keep the pointer to mpm
  string equation;             // formula
  double value;                // current value
  bool constant;               // is the variables constant?
};


Var powv(int, Var);            // power function

#endif
