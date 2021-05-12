#include <iostream>
#include "temperature.h"

using namespace std;


Temperature::Temperature(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new Temperature with ID: " << args[0] << endl;
  id = args[0];
  style = args[1];
}

Temperature::~Temperature()
{
}


void Temperature::init()
{
}

void Temperature::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In Temperature::options()" << endl;
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

