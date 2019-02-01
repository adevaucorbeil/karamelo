#include <iostream>
#include "domain.h"
#include "region_block.h"

using namespace std;

RegBlock::RegBlock(MPM *mpm, string *args) : Region(mpm, args)
{
  cout << "Initiate RegBlock" << endl; 
}


RegBlock::~RegBlock()
{

}
