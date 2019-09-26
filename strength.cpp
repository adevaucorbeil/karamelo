#include <iostream>
#include "strength.h"
#include "error.h"

using namespace std;


Strength::Strength(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new Strength with ID: " << args[0] << endl;
  id = args[0];
}

Strength::~Strength()
{
}


void Strength::init()
{
}

void Strength::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In Strength::options()" << endl;
  if (args->end() < it) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

