#include "method.h"

using namespace std;


Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_total_lagrangian = false;
}

Method::~Method()
{
}
