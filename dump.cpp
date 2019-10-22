#include <iostream>
#include <algorithm>
#include "dump.h"

using namespace std;


Dump::Dump(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new dump with ID: " << args[0] << endl;
  id = args[0];

  if (args[1].compare("all")!=0) {
    cout << "Error: groups are not yet supported for dumps. Should use all." << endl;
    exit(1);
  }

  style = args[2];
  
  filename = args[4];
}

Dump::~Dump()
{
}
