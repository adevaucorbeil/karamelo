#include <iostream>
#include "compute.h"

using namespace std;


Compute::Compute(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new compute with ID: " << args[0] << endl;
  id = args[0];
}
