/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_VARIABLE_H
#define MPM_VARIABLE_H

#include <iostream>
#include <string>
#include <vector>

#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;
class Variable {
public:
  map<string, class Variable> *known_var;

  Variable() {};
  Variable(map<string, class Variable> *, double, bool c = true);
  Variable(map<string, class Variable> *, string, double, bool c = false);
  ~Variable() {};

  void evaluate();
  double result();
  double result() const;
  string str() const;
  string eq() const;
  bool is_constant() const;
  Variable operator+(const Variable&);
  Variable operator-(const Variable&);
  Variable operator-();
  Variable operator*(const Variable&);
  Variable operator/(const Variable&);
  Variable operator^(const Variable&);

protected:
  string equation;
  double value;
  bool constant;

};

double precedence(char);
double applyOp(double, double, char);
Variable applyOp(Variable, Variable, char);
bool is_operator(char);
bool is_math_char(char);
Variable evaluate_function(map<string, Variable> *, string, string);
string remove_whitespace(string);
Variable parse(map<string, Variable> *, string);
Variable powv(int, Variable);
Variable powv(Variable, Variable);
Variable expv(Variable);

Variable operator*(int, Variable);
Variable operator*(Variable, int);

#endif
