#include <iostream>
#include "group.h"
#include "var.h"
#include "centre_of_mass.h"

using namespace std;

CentreOfMass::CentreOfMass(MPM *mpm) : Pointers(mpm) {}

Var CentreOfMass::command(vector<string> args) {
  // cout << "In CentreOfMass::command()" << endl;

  if (args.size() < 2) {
    cout << "Illegal run command" << endl;
    exit(1);
  }

  igroup = group->find(args[0]);

  if (igroup == -1) {
    cout << "Error: could not find group named: " << args[0] << endl;
  }

  if (args[1].compare("x") == 0) dir = 0;
  else if (args[1].compare("y") == 0) dir = 1;
  else if (args[1].compare("z") == 0) dir = 2;
  else {
    cout << "Error: directions should be either x,y or z: " << args[1] << " not understood." << endl;
    exit(1);
  }

  double com = group->xcm(igroup, dir);
  string eq = "xcm(" + args[0] + "," + args[1] + ")";
  return Var(eq, com);
}
