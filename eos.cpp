#include <iostream>
#include "eos.h"
#include "error.h"

using namespace std;


EOS::EOS(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new EOS with ID: " << args[0] << endl;
  id = args[0];
}

EOS::~EOS()
{
}


void EOS::init()
{
}

void EOS::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In EOS::options()" << endl;
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

