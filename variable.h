/* -*- c++ -*- ----------------------------------------------------------*/

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Variable;

map<string, Variable> variables; // global variables

class Variable {
public:

  
  Variable() {};
  Variable(double, bool c = true);
  Variable(string, double, bool c = false);
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

double precedence(char op);
double applyOp(double a, double b, char op);
Variable applyOp(Variable a, Variable b, char op);
bool is_operator(char op);
bool is_math_char(char op);
Variable evaluate_function(string func, string arg);
string remove_whitespace(string str);
Variable parse(string str);
Variable powv(int, Variable);
Variable powv(Variable, Variable);
Variable expv(Variable);

Variable operator*(int, Variable);
Variable operator*(Variable, int);
