#include <iostream>
#include "log.h"

using namespace std;

Log::Log(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  cout << "In Log::Log" << endl;

  style = args[0];

  if (style.compare("default")){
    
  } else {
    cout << "Unknown log style " << style << endl;
    exit(1);
  }
}

void Log::write(){
  cout << "Writing output ... function to write" << endl;
}
