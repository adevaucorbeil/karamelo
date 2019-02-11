#include <iostream>
#include <fstream>
#include "log.h"
#include "mpmtype.h"
#include "update.h"
#include <sstream>

using namespace std;

enum{INT,FLOAT,BIGINT};

Log::Log(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  ivalue = 0;
  dvalue = 0;
  bivalue = 0;

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
  stringstream soutput;
  string output = "";

  for (int i=0; i<field.size(); i++){
    (this->*field[i].vfunc)(); // Compute the output field

    //if (field[i].typeflag==INT) cout << to_string(ivalue) << "\t";
    //if (field[i].typeflag==FLOAT) cout << to_string(dvalue) << "\t";
    //if (field[i].typeflag==BIGINT) cout << to_string(bivalue) << "\t";
    soutput << ivalue << "\t";
    soutput << dvalue << "\t";
    soutput << bivalue << "\t";
  }

  soutput << "\n";
  cout << soutput.str();

  if (wlogfile->is_open()) {
    (*wlogfile) << soutput.str();
  }
}

void Log::parse_keywords(vector<string> keyword)
{
  for (int i=0; i<keyword.size();i++){
    if (keyword[i].compare("step")==0)      addfield("Step", &Log::compute_step, BIGINT);
    else if (keyword[i].compare("dt")==0)   addfield("dt", &Log::compute_dt, FLOAT);
    else if (keyword[i].compare("time")==0) addfield("Time", &Log::compute_time, FLOAT);
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
  Field f = {};
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
