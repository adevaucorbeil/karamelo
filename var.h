/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_VAR_H
#define MPM_VAR_H

#include "input.h"
#include <vector>

class Var : protected Input {
 public:
  // functions
  Var(class MPM *) : Input(mpm, 0, NULL) {};
  Var() : Input(NULL, 0, NULL) {};
  Var(class MPM*, double, bool c = true);
  Var(class MPM*, string, double, bool c = false);
  ~Var() {};

  void evaluate();
  double result();
  double result() const {return value;};
  string str() const;
  string eq() const;
  bool is_constant() const {return constant;};
  Var operator+(const Var&);
  Var operator-(const Var&);
  Var operator-();
  Var operator*(const Var&);
  Var operator/(const Var&);
  Var operator^(const Var&);
  Var operator=(const Var&);

protected:
  string equation;
  double value;
  bool constant;
};

#endif
