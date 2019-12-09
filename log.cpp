#include <iostream>
#include <fstream>
#include "log.h"
#include "mpmtype.h"
#include "update.h"
#include "input.h"
#include "var.h"
#include <sstream>
#include <stdexcept>

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
    (this->*field[i].vfunc)(this->field[i].name); // Compute the output field

    if (field[i].typeflag==INT) soutput << ivalue << "\t";
    if (field[i].typeflag==FLOAT) soutput << dvalue << "\t";
    if (field[i].typeflag==BIGINT) soutput << bivalue << "\t";
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
      try {
	(*input->vars).at(keyword[i]);
	addfield(keyword[i], &Log::compute_var, FLOAT);
      }
      catch (const std::out_of_range& oor) {
	cout << "Error: unknown log keyword " << keyword[i] << endl;
	std::cerr << "Out of Range error: " << oor.what() << '\n';
	exit(1);
      }
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

void Log::clearfields()
{
  field.clear();
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

void Log::compute_step(string name)
{
  bivalue = update->ntimestep;
}

void Log::compute_dt(string name)
{
  dvalue = update->dt;
}

void Log::compute_time(string name)
{
  dvalue = update->atime;// + (update->ntimestep-update->atimestep)*update->dt;
}

void Log::compute_var(string name)
{
  dvalue = (*input->vars)[name].result(mpm);
}

void Log::modify(vector<string> args)
{
  if (args.size() < 1) {
    cout << "Erro: too few arguments given to log_modify command." << endl;
    exit(1);
  }

  style = args[0];

  vector<string> keyword;

  if (style.compare("default") == 0){
    keyword.push_back("step");
    keyword.push_back("dt");
    keyword.push_back("time");
  } else if (style.compare("custom") == 0) {
    
    if (args.size() < 2) {
      cout << "Error: too few arguments given to log_modify(custom,...) command." << endl;
      exit(1);
    }

    for (int i=1; i<args.size(); i++) keyword.push_back(args[i]);

  } else {
    cout << "Unknown log style " << style << endl;
    exit(1);
  }

  clearfields();
  parse_keywords(keyword);
}
