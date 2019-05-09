#include "mpm.h"
#include "input.h"
#include "var.h"
#include "math.h"

using namespace std;

Var::Var(double v)
{
  equation = to_string(v);
  value = v;
  constant = true;
}

Var::Var(string eq, double v, bool c)
{
  equation = eq;
  value = v;
  constant = c;
}

Var::~Var()
{
}

void Var::evaluate(MPM * mpm)
{
  if (constant) return;
  if (mpm) value = mpm->input->parsev(equation).value;
  if (constant) equation = to_string(value);
}


double Var::result(MPM * mpm)
{
  if (constant) return value;
  else {
    evaluate(mpm);
    return value;
  }
}

string Var::str() const
{
  if (equation != "") return equation;
  else return to_string(value);
}

void Var::make_constant(MPM * mpm)
{
  if (!constant) {
    if (mpm) value = mpm->input->parsev(equation).value;
    equation = to_string(value);
    constant = true;
  }
}

// Var Var::operator=(const Var& right){
//   return right;
//}


Var Var::operator+(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value + right.value);
    return result;
  } else {
    Var result("(" + this->str() + "+" + right.str() + ")", this->value + right.value, this->constant && right.constant);
    return result;
  }
}

Var Var::operator-(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value - right.value);
    return result;
  } else {
    Var result("(" + this->str() + "-" + right.str() + ")", this->value - right.value, this->constant && right.constant);
    return result;
  }
}

Var Var::operator-()
{
  if (this->constant) {
    Var result(-this->value);
    return result;
  } else {
    Var result("(-" + this->str() + ")", -this->value, this->constant);
    return result;
  }
}

Var Var::operator*(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value * right.value);
    return result;
  } else {
    Var result("(" + this->str() + "*" + right.str() + ")", this->value * right.value, false);
    return result;
  }
}


Var Var::operator/(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value / right.value);
    return result;
  } else {
    Var result("(" + this->str() + "/" + right.str() + ")", this->value / right.value, false);
    return result;
  }
}


Var Var::operator^(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(pow(this->value, right.value));
    return result;
  } else {
    Var result("(" + this->str() + "^" + right.str() + ")", pow(this->value, right.value), false);
    return result;
  }
}

Var powv(int base, Var p){
  if (p.is_constant()) {
    Var result(pow(base, p.result()));
    return result;
  } else {
    Var result("pow(" + to_string(base) + "," + p.str() + ")", pow(base, p.result()), p.is_constant());
    return result;
  }
}

Var expv(Var x){
  if (x.is_constant()) {
    Var result(exp(x.result()));
    return result;
  } else {
    Var result("exp(" + x.str() + ")", exp(x.result()), x.is_constant());
    return result;
  }
}


Var sqrtv(Var x){
  if (x.is_constant()) {
    Var result(sqrt(x.result()));
    return result;
  } else {
    Var result("sqrt(" + x.str() + ")", sqrt(x.result()), x.is_constant());
    return result;
  }
}


Var cosv(Var x){
  if (x.is_constant()) {
    Var result(cos(x.result()));
    return result;
  } else {
    Var result("cos(" + x.str() + ")", cos(x.result()), x.is_constant());
    return result;
  }
}


Var sinv(Var x){
  if (x.is_constant()) {
    Var result(sin(x.result()));
    return result;
  } else {
    Var result("sin(" + x.str() + ")", sin(x.result()), x.is_constant());
    return result;
  }
}

Var logv(Var x){
  if (x.is_constant()) {
    Var result(log(x.result()));
    return result;
  } else {
    Var result("log(" + x.str() + ")", log(x.result()), x.is_constant());
    return result;
  }
}
