#include "mpm.h"
#include "input.h"
#include "var.h"
#include "math.h"

using namespace std;

Var::Var() : Input(NULL, 0, NULL)
{
  mpm_ptr = NULL;
}

Var::Var(MPM *mpm) : Input(mpm, 0, NULL)
{
  mpm_ptr = mpm;
}


Var::Var(MPM *mpm, double v, bool c) : Input(mpm, 0, NULL)
{
  mpm_ptr = mpm;
  equation = to_string(v);
  value = v;
  constant = c;
}

Var::Var(MPM *mpm, string eq, double v, bool c) : Input(mpm, 0, NULL)
{
  mpm_ptr = mpm;
  equation = eq;
  value = v;
  constant = c;
}

void Var::evaluate()
{
  if (constant) return;
  value = input->parsev(equation).value;
  if (constant) equation = to_string(value);
}


double Var::result()
{
  if (constant) return value;
  else {
    evaluate();
    return value;
  }
}

string Var::str() const
{
  if (equation != "") return equation;
  else return to_string(value);
}

MPM* Var::mpm() const
{
  return mpm_ptr;
}

Var Var::operator=(const Var& right){
  Var result(right.mpm());
  return result;
}


Var Var::operator+(const Var& right)
{
  if (this->mpm() != right.mpm()) {
    cout << "Error: this->mpm() != right.mpm(), I don't know how to deal with this" << endl;
    exit(1);
  }
  if (this->constant && right.constant) {
    Var result(this->mpm(), this->value + right.value);
    return result;
  } else {
    Var result(this->mpm(), "(" + this->str() + "+" + right.str() + ")", this->value + right.value, this->constant && right.constant);
    return result;
  }
}

Var Var::operator-(const Var& right)
{
  if (this->mpm() != right.mpm()) {
    cout << "Error: this->mpm() != right.mpm(), I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Var result(this->mpm(), this->value - right.value);
    return result;
  } else {
    Var result(this->mpm(), "(" + this->str() + "-" + right.str() + ")", this->value - right.value, this->constant && right.constant);
    return result;
  }
}

Var Var::operator-()
{
  if (this->constant) {
    Var result(this->mpm(), -this->value);
    return result;
  } else {
    Var result(this->mpm(), "(-" + this->str() + ")", -this->value, this->constant);
    return result;
  }
}

Var Var::operator*(const Var& right)
{
  if (this->mpm() != right.mpm()) {
    cout << "Error: this->mpm() != right.mpm(), I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Var result(this->mpm(), this->value * right.value);
    return result;
  } else {
    Var result(this->mpm(), "(" + this->str() + "*" + right.str() + ")", this->value * right.value, this->constant && right.constant);
    return result;
  }
}


Var Var::operator/(const Var& right)
{
  if (this->mpm() != right.mpm()) {
    cout << "Error: this->mpm() != right.mpm(), I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Var result(this->mpm(), this->value / right.value);
    return result;
  } else {
    Var result(this->mpm(), "(" + this->str() + "/" + right.str() + ")", this->value / right.value, this->constant && right.constant);
    return result;
  }
}


Var Var::operator^(const Var& right)
{
  if (this->mpm() != right.mpm()) {
    cout << "Error: this->mpm() != right.mpm(), I don't know how to deal with this" << endl;
    exit(1);
  }

  if (this->constant && right.constant) {
    Var result(this->mpm(), pow(this->value, right.value));
    return result;
  } else {
    Var result(this->mpm(), "(" + this->str() + "^" + right.str() + ")", pow(this->value, right.value), this->constant && right.constant);
    return result;
  }
}

Var powv(int base, Var p){
  if (p.is_constant()) {
    Var result(p.mpm(), pow(base, p.result()));
    return result;
  } else {
    Var result(p.mpm(), "pow(" + to_string(base) + "," + p.str() + ")", pow(base, p.result()), p.is_constant());
    return result;
  }
}

Var expv(Var x){
  if (x.is_constant()) {
    Var result(x.mpm(), exp(x.result()));
    return result;
  } else {
    Var result(x.mpm(), "exp(" + x.str() + ")", exp(x.result()), x.is_constant());
    return result;
  }
}
