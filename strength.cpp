#include <iostream>
#include "strength.h"

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

