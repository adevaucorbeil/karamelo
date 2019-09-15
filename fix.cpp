#include <iostream>
#include "fix.h"
#include "group.h"
#include "input.h"
#include "var.h"

using namespace std;


Fix::Fix(MPM *mpm, vector<string> args_) :
  Pointers(mpm)
{
  mask = 0;
  args = args_;
  cout << "Creating new fix with ID: " << args[0] << endl;
  id = args[0];

  igroup = group->find(args[2]);
  if (igroup == -1) {
    cout << "Could not find group ID " << args[2] << endl;
  }
  groupbit = group->bitmask[igroup];
  (*input->vars)[id+"_x"]=Var(id+"_x", 0);
  (*input->vars)[id+"_y"]=Var(id+"_y", 0);
  (*input->vars)[id+"_z"]=Var(id+"_z", 0);
  (*input->vars)[id+"_s"]=Var(id+"_s", 0);

}
