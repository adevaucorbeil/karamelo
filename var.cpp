#include "mpm.h"
#include "input.h"
#include "var.h"

using namespace std;

Var::Var(MPM *mpm, double v, bool c) : Input(mpm, 0, NULL)
{
  equation = to_string(v);
  value = v;
  constant = c;
}

Var::Var(MPM *mpm, string eq, double v, bool c) : Input(mpm, 0, NULL)
{
  equation = eq;
  value = v;
  constant = c;
}

Var Var::operator=(const Var& right){
  Var result(right.mpm);
  return result;
}
