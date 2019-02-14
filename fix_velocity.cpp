#include <iostream>
#include <vector>
#include "fix_velocity.h"

using namespace std;


FixVelocity::FixVelocity(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  cout << "Creating new fix FixVelocity with ID: " << args[0] << endl;
  id = args[0];
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
