#include <iostream>
#include "fix.h"
#include "group.h"

using namespace std;


Fix::Fix(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new fix with ID: " << args[0] << endl;
  id = args[0];

  igroup = group->find(args[2]);
  if (igroup == -1) {
    cout << "Could not find group ID " << args[2] << endl;
  }
  groupbit = group->bitmask[igroup];

}
