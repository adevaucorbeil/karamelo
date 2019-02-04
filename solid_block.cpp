#include <iostream>
#include "domain.h"
#include "solid_block.h"
#include "input.h"

using namespace std;


SolBlock::SolBlock(MPM *mpm, vector<string> args) : Solid(mpm, args)
{
  cout << "Initiate SolBlock" << endl;

  if (args.size() < domain->dimension + 3) {
    cout << "Error: solid command not enough arguments" << endl;
    exit(1);
  }
  options(&args, args.begin()+domain->dimension + 3);

  cout << "Solid delimitated by region ID: " << args[2] << endl;

  // Look for region ID:
  int iregion = domain->find_region(args[2]);
  if (iregion == -1) {
    cout << "Error: region ID " << args[2] << " not does not exist" << endl;
    exit(1);
  }
  
}


SolBlock::~SolBlock()
{

}
