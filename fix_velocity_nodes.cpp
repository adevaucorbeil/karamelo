#include <iostream>
#include <vector>
#include "fix_velocity_nodes.h"
#include "input.h"
#include "group.h"

using namespace std;
using namespace FixConst;

FixVelocityNodes::FixVelocityNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_velocity_nodes: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    cout << "fix_velocity_nodes needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixVelocityNodes with ID: " << args[0] << endl;
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

FixVelocityNodes::~FixVelocityNodes()
{
}

void FixVelocityNodes::init()
{
}

void FixVelocityNodes::setup()
{
}

void FixVelocityNodes::setmask() {
  mask = 0;
  // mask |= POST_UPDATE_GRID_STATE;
  // mask |= POST_VELOCITIES_TO_GRID;
}


void FixVelocityNodes::post_update_grid_state() {
  cout << "In FixVelocityNodes::post_update_grid_state()" << endl;

  // Go through all the particles in the group and set their velocities 
}

void FixVelocityNodes::post_velocities_to_grid() {
  cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;
}
