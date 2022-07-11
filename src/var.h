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
#include <string>

using namespace std;

class Var{
 public:
  // functions
  Var() {};
  Var(float);
  Var(string, float, bool c = false);

  void evaluate(class MPM * = nullptr);
  float result(class MPM *mpm, bool only_position_changed = false);
  float result() const {return value;};
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
  operator float() {return value;};

  void write_to_restart(ofstream *);
  void read_from_restart(ifstream *);

protected:
  string equation;             // formula
  float value;                // current value
  bool constant;               // is the variables constant?
  bool position_independent = false;
};


Var powv(int, Var);            // power function
Var expv(Var);                 // exponential function
Var sqrtv(Var);                // square root function
Var cosv(Var);                 // cosine function
Var sinv(Var);                 // sine function
Var tanv(Var);                 // tangent function
Var logv(Var);                 // natural logarithmic function
Var atan2v(Var, Var);          // arc tangent function of two variables
Var ifv(Var, Var, Var);        // if statement

#endif
