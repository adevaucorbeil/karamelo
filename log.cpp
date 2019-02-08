#include <iostream>
#include <fstream>
#include "log.h"
#include "mpmtype.h"
#include "update.h"

using namespace std;

enum{INT,FLOAT,BIGINT};

  Log::Log(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  cout << "In Log::Log" << endl;

  style = args[0];

  vector<string> keyword;

  if (style.compare("default") == 0){
    keyword.push_back("step");
    keyword.push_back("dt");
    keyword.push_back("time");
  } else {
    cout << "Unknown log style " << style << endl;
    exit(1);
  }
  parse_keywords(keyword);
}

void Log::write()
{
  cout << "Writing output ... function to write" << endl;
}

void Log::parse_keywords(vector<string> keyword)
{
  for (int i=0; i<keyword.size();i++){
    if (keyword[i].compare("step")==0)      addfield("Step", &Log::compute_step, BIGINT);
    else if (keyword[i].compare("dt")==0)   addfield("dt", &Log::compute_dt, BIGINT);
    else if (keyword[i].compare("time")==0) addfield("Time", &Log::compute_time, BIGINT);
    else {
      cout << "Error: unknown log keyword " << keyword[i] << endl;
      exit(1);
    }
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Log::addfield(string key, FnPtr func, int typeflag)
{
  Field f;
  f.name = key;
  f.vfunc = func;
  f.typeflag = typeflag;
  field.push_back(f);
}

void Log::init()
{
}

void Log::header()
{
  string output = "";
  for (int i=0; i<field.size(); i++){
    output += field[i].name + "\t";
  }
  output += "\n";
  cout << output;
  if (wlogfile->is_open()) {
    (*wlogfile) << output;
  }
}

void Log::compute_step()
{
  bivalue = update->ntimestep;
}

void Log::compute_dt()
{
  dvalue = update->dt;
}

void Log::compute_time()
{
  dvalue = update->atime + (update->ntimestep-update->atimestep)*update->dt;
}
