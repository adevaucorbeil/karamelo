#include <iostream>
#include "output.h"
#include "dump_particle.h"

using namespace std;


DumpParticle::DumpParticle(MPM *mpm, vector<string> args) : Dump(mpm, args)
{
  cout << "In DumpParticle::DumpParticle()" << endl;
}

DumpParticle::~DumpParticle()
{
}
