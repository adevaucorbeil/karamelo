#include "mpm.h"
#include "solid.h"
#include <vector>

using namespace std;


Solid::Solid(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new solid with ID: " << args[0] << endl;
  id = args[0];
}

Solid::~Solid()
{
}


void Solid::init()
{
}

void Solid::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In solid::options()" << endl;
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
