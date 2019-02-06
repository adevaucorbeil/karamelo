#include <iostream>
#include <vector>
#include "run.h"
#include "domain.h"
#include "output.h"
#include "input.h"


/* ---------------------------------------------------------------------- */

Run::Run(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

void Run::command(vector<string> args)
{
  cout << "In Run::command()" << endl;
  if (args.size() < 1) {
    cout << "Illegal run command" << endl;
    exit(1);
  }
}
