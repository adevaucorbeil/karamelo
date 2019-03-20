#include <iostream>
#include "group.h"
#include "var.h"
#include "internal_force.h"

using namespace std;

IntForce::IntForce(MPM *mpm) : Pointers(mpm) {}

Var IntForce::command(vector<string> args) {
  // cout << "In IntForce::command()" << endl;

  if (args.size() < 2) {
    cout << "Illegal run command" << endl;
    exit(1);
  }

  igroup = group->find(args[0]);

  if (igroup == -1) {
    cout << "Error: could not find group named: " << args[0] << endl;
    exit(1);
  }

  if (args[1].compare("x") == 0) dir = 0;
  else if (args[1].compare("y") == 0) dir = 1;
  else if (args[1].compare("z") == 0) dir = 2;
  else {
    cout << "Error: directions should be either x,y or z: " << args[1] << " not understood." << endl;
    exit(1);
  }

  double iforce = group->internal_force(igroup, dir);
  string eq = "internal_force(" + args[0] + "," + args[1] + ")";
  return Var(eq, iforce);
}
