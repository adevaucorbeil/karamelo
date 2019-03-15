#include <iostream>
#include <vector>
#include "fix_velocity.h"
#include "input.h"

using namespace std;


FixVelocity::FixVelocity(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_velocity: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }
  cout << "Creating new fix FixVelocity with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  if (args[3].compare("NULL") != 0) {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (args[4].compare("NULL") != 0) {
    yvalue = input->parsev(args[4]);
    yset = true;
  }

  if (args[5].compare("NULL") != 0) {
    yvalue = input->parsev(args[5]);
    zset = true;
  }
}

FixVelocity::~FixVelocity()
{
}

void FixVelocity::init()
{
}

void FixVelocity::setup()
{
}
