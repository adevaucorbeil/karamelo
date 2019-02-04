#include <iostream>
#include "domain.h"
#include "solid_block.h"
#include "input.h"

using namespace std;


SolBlock::SolBlock(MPM *mpm, vector<string> args) : Solid(mpm, args)
{
  cout << "Initiate SolBlock" << endl;

  // if (args.size()<8) {
  //   cout << "Error: solid command not enough arguments" << endl;
  //   exit(1);
  // }
  // options(&args, args.begin()+8);

}


SolBlock::~SolBlock()
{

}
