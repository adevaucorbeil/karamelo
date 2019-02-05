#include <iostream>
#include "mat.h"

using namespace std;


Mat::Mat(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new mat with ID: " << args[0] << endl;
  id = args[0];
}

Mat::~Mat()
{
}


void Mat::init()
{
}

void Mat::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In mat::options()" << endl;
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

