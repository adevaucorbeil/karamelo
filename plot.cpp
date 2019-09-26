#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <matplotlibcpp.h>
#include "plot.h"
#include "mpmtype.h"
#include "update.h"
#include "input.h"
#include "var.h"
#include "error.h"

namespace plt = matplotlibcpp;
using namespace std;

enum{INT,FLOAT,BIGINT};

Plot::Plot(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  id = args[0];

  ivalue = 0;
  dvalue = 0;
  bivalue = 0;

  cout << "In Plot::Plot" << endl;

  style = args[0];

  vector<string> keyword;

  keyword.push_back(args[2]); // x
  keyword.push_back(args[3]); // y
  parse_keywords(keyword);
}

void Plot::save()
{
  for (int i=0; i<field.size(); i++){
    (this->*field[i].vfunc)(this->field[i].name, this->field[i].x_or_y); // Compute the output field
  }
}

void Plot::parse_keywords(vector<string> keyword)
{
  int x_or_y = 0;
  for (int i=0; i<keyword.size();i++){

    if (i==0) x_or_y = 0;
    else if (i==1) x_or_y = 1;
    else {
      error->all(FLERR, "Error: too many variables to plot.\n");
    }

    if (keyword[i].compare("step")==0)      addfield("Step", &Plot::compute_step, BIGINT, x_or_y);
    else if (keyword[i].compare("dt")==0)   addfield("dt", &Plot::compute_dt, FLOAT, x_or_y);
    else if (keyword[i].compare("time")==0) addfield("Time", &Plot::compute_time, FLOAT, x_or_y);
    else {
      try {
	(*input->vars).at(keyword[i]);
	addfield(keyword[i], &Plot::compute_var, FLOAT, x_or_y);
      }
      catch (const std::out_of_range& oor) {
	error->all(FLERR, "Error: unknown plot keyword " + keyword[i] + ".\n");
	// std::cerr << "Out of Range error: " << oor.what() << '\n';
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Plot::addfield(string key, FnPtrPlt func, int typeflag, int x_or_y)
{
  FieldPlt f = {};
  f.name = key;
  f.vfunc = func;
  f.typeflag = typeflag;
  f.x_or_y = x_or_y;
  field.push_back(f);
}

void Plot::clearfields()
{
  field.clear();
}

void Plot::init()
{
}

void Plot::compute_step(string name, int x_or_y)
{
  if (x_or_y == 0) {
    x.push_back(update->ntimestep);
  } else if (x_or_y == 1) {
    y.push_back(update->ntimestep);
  } else {
    error->all(FLERR, "Error: in Plot::compute_step, x_or_y == " + to_string(x_or_y) + ". Expected 0 or 1!\n");
  }
}

void Plot::compute_dt(string name, int x_or_y)
{
  if (x_or_y == 0) {
    x.push_back(update->dt);
  } else if (x_or_y == 1) {
    y.push_back(update->dt);
  } else {
    error->all(FLERR, "Error: in Plot::compute_dt, x_or_y == " + to_string(x_or_y) + ". Expected 0 or 1!\n");
  }
}

void Plot::compute_time(string name, int x_or_y)
{
  if (x_or_y == 0) {
    x.push_back(update->atime);
  } else if (x_or_y == 1) {
    y.push_back(update->atime);
  } else {
    error->all(FLERR, "Error: in Plot::compute_time, x_or_y == " + to_string(x_or_y) + ". Expected 0 or 1!\n");
  }
}

void Plot::compute_var(string name, int x_or_y)
{
  if (x_or_y == 0) {
    x.push_back((*input->vars)[name].result(mpm));
  } else if (x_or_y == 1) {
    y.push_back((*input->vars)[name].result(mpm));
  } else {
    error->all(FLERR, "Error: in Plot::compute_var, x_or_y == " + to_string(x_or_y) + ". Expected 0 or 1!\n");
  }
}

void Plot::modify(vector<string> args)
{
  if (args.size() < 1) {
    error->all(FLERR, "Error: too few arguments given to plot_modify command.\n");
  }

  style = args[0];

  vector<string> keyword;

  keyword.push_back(args[2]); // x
  keyword.push_back(args[3]); // y
  clearfields();
  parse_keywords(keyword);
}

