/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#ifndef MPM_VAR_H
#define MPM_VAR_H

#include <vector>

class Var{
 public:
  // functions
  Var() {};
  Var(double);
  Var(string, double, bool c = false);

  void evaluate(class MPM * = nullptr);
  double result(class MPM *);
  double result() const {return value;};
  string str() const;
  string eq() const {return equation;};
  bool is_constant() const {return constant;};
  void make_constant(class MPM *);
  Var operator+(const Var&);
  Var operator-(const Var&);
  Var operator-();
  Var operator*(const Var&);
  Var operator/(const Var&);
  Var operator^(const Var&);
  Var operator>(const Var&);
  Var operator>=(const Var&);
  Var operator<(const Var&);
  Var operator<=(const Var&);
  Var operator==(const Var&);
  Var operator!=(const Var&);
  Var operator!();
  //Var operator=(const Var&);
  operator int() {return (int) value;};
  operator double() {return value;};

  void write_to_restart(ofstream *);
  void read_from_restart(ifstream *);

protected:
  string equation;             // formula
  double value;                // current value
  bool constant;               // is the variables constant?
};


Var powv(int, Var);            // power function
Var expv(Var);                 // exponential function
Var sqrtv(Var);                // square root function
Var cosv(Var);                 // cosine function
Var sinv(Var);                 // sine function
Var tanv(Var);                 // tangent function
Var logv(Var);                 // natural logarithmic function
Var atan2v(Var, Var);          // arc tangent function of two variables

#endif
