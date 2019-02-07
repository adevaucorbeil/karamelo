#include <iostream>
#include "region.h"

using namespace std;


Region::Region(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new region with ID: " << args[0] << endl;
  id = args[0];
}

Region::~Region()
{
}


void Region::init()
{
}

void Region::options(vector<string> *args, vector<string>::iterator it)
{
  interior = 1; // Interior by default

  cout << "In region::options()" << endl;
  if (args->end() < it) {
    cout << "Error: not enough arguments" << endl;
    exit(1);
  }
  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

int Region::match(double x, double y, double z)
{
  return !(inside(x,y,z) ^ interior);
}
