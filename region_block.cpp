#include <iostream>
#include "domain.h"
#include "region_block.h"

using namespace std;

RegBlock::RegBlock(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegBlock" << endl; 
}


RegBlock::~RegBlock()
{

}
