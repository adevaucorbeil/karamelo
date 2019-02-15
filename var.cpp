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

Var Var::operator=(const Var& right){
  Var result(right.mpm);
  return result;
}
