#include <iostream>
#include "fix.h"

using namespace std;


Fix::Fix(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new fix with ID: " << args[0] << endl;
  id = args[0];
}
