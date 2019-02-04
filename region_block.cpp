#include <iostream>
#include "domain.h"
#include "region_block.h"
#include "input.h"

using namespace std;

#define BIG 1.0e20


RegBlock::RegBlock(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegBlock" << endl;

  if (args.size()<8) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }
  options(&args, args.begin()+8);

  xlo = input->parse(args[2]);
  xhi = input->parse(args[3]);
  cout << "xlo xhi = " << xlo << "\t" << xhi << endl;
  ylo = input->parse(args[4]);
  yhi = input->parse(args[5]);
  cout << "ylo yhi = " << ylo << "\t" << yhi << endl;
  zlo = input->parse(args[6]);
  zhi = input->parse(args[7]);
  cout << "zlo zhi = " << zlo << "\t" << zhi << endl;
}


RegBlock::~RegBlock()
{

}
