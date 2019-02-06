#include "musl.h"
#include <iostream>
#include <vector>

using namespace std;

MUSL::MUSL(MPM *mpm, vector<string> args) : Scheme(mpm) {
  cout << "In MUSL::MUSL()" << endl;
}

void MUSL::run(int nsteps){
  cout << "In MUSL::run" << endl;
}

