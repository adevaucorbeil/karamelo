#include "tlmpm.h"
#include <iostream>
#include <vector>

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;
}

