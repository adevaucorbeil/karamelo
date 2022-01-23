/* ----------------------------------------------------------------------
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

#include <mpm.h>
#include <input.h>
#include <var.h>
#include <math.h>
#include <sstream>

using namespace std;

//#define PRECISION 15

std::string to_string_with_precision(const double a_value){
    std::string s(30, '\0');
    auto written = std::snprintf(&s[0], s.size(), "%.15f", a_value);
    s.resize(written);
    return s;
}

Var::Var(double v)
{
  equation = to_string_with_precision(v);
  value = v;
  constant = true;
}

Var::Var(string eq, double v, bool c)
{
  equation = eq;
  value = v;
  constant = c;
}

void Var::evaluate(MPM * mpm)
{
  if (constant) return;
  if (mpm) value = mpm->input->parsev(equation).value;
  if (constant) equation = to_string_with_precision(value);
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
  else return to_string_with_precision(value);
}

void Var::make_constant(MPM * mpm)
{
  if (!constant) {
    if (mpm) value = mpm->input->parsev(equation).value;
    equation = to_string_with_precision(value);
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
    return Var(this->value * right.value);
    //return result;
  } else {
    return Var("(" + this->str() + "*" + right.str() + ")", this->value * right.value, false);
    //return result;
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

Var Var::operator>(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value > right.value);
    return result;
  } else {
    Var result("(" + this->str() + ">" + right.str() + ")", this->value > right.value, false);
    return result;
  }
}

Var Var::operator>=(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value >= right.value);
    return result;
  } else {
    Var result("(" + this->str() + ">=" + right.str() + ")", this->value >= right.value, false);
    return result;
  }
}

Var Var::operator<(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value < right.value);
    return result;
  } else {
    Var result("(" + this->str() + "<" + right.str() + ")", this->value < right.value, false);
    return result;
  }
}

Var Var::operator<=(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value <= right.value);
    return result;
  } else {
    Var result("(" + this->str() + "<=" + right.str() + ")", this->value <= right.value, false);
    return result;
  }
}

Var Var::operator==(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value == right.value);
    return result;
  } else {
    Var result("(" + this->str() + "==" + right.str() + ")", this->value == right.value, false);
    return result;
  }
}

Var Var::operator!=(const Var& right)
{
  if (this->constant && right.constant) {
    Var result(this->value != right.value);
    return result;
  } else {
    Var result("(" + this->str() + "!=" + right.str() + ")", this->value != right.value, false);
    return result;
  }
}

Var Var::operator!()
{
  if (this->constant) {
    Var result(!this->value);
    return result;
  } else {
    Var result("(!" + this->str() + ")", !this->value, this->constant);
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

Var tanv(Var x){
  if (x.is_constant()) {
    Var result(tan(x.result()));
    return result;
  } else {
    Var result("tan(" + x.str() + ")", tan(x.result()), x.is_constant());
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

Var atan2v(Var x, Var y){
  if (x.is_constant() && y.is_constant()) {
    Var result(atan2(x.result(), y.result()));
    return result;
  } else {
    Var result("atan2(" + x.str() + ", " + y.str() + ")", atan2(x.result(), y.result()), false);
    return result;
  }
}

Var ifv(Var condition, Var x, Var y) {
  if (condition.is_constant()) {
    if (condition.result())
      return x;
    else
      return y;
  } else {
    double value;
    if (condition.result())
      value = x.result();
    else
      value = y.result();
    Var result("if(" + condition.str() + "," + x.str() + "," + y.str() + ")",
               value, false);
    return result;
  }
}

void Var::write_to_restart(ofstream *of) {
  size_t N = equation.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  of->write(reinterpret_cast<const char *>(equation.c_str()), N);
  of->write(reinterpret_cast<const char *>(&value), sizeof(double));
  of->write(reinterpret_cast<const char *>(&constant), sizeof(bool));

}


void Var::read_from_restart(ifstream *ifr) {
  size_t N = 0;

  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  equation.resize(N);

  ifr->read(reinterpret_cast<char *>(&equation[0]), N);
  ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&constant), sizeof(bool));
}
