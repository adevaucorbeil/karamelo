#include "method.h"

using namespace std;


Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
}

Method::~Method()
{
}
