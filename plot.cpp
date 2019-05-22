#include <iostream>
#include <fstream>
#include "plot.h"
#include "mpmtype.h"
#include "update.h"
#include "input.h"
#include "var.h"
#include <sstream>
#include <stdexcept>
#include <matplotlibcpp.h>

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
      cout << "Error: too many variables to plot." << endl;
      exit(1);
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
	cout << "Error: unknown plot keyword " << keyword[i] << endl;
	std::cerr << "Out of Range error: " << oor.what() << '\n';
	exit(1);
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
    xbi.push_back(update->ntimestep);
  } else if (x_or_y == 1) {
    ybi.push_back(update->ntimestep);
  } else {
    cout << "Error: in Plot::compute_step, x_or_y == " << x_or_y << ". Expected 0 or 1!" << endl;
    exit(1);
  }
}

void Plot::compute_dt(string name, int x_or_y)
{
  if (x_or_y == 0) {
    xd.push_back(update->dt);
  } else if (x_or_y == 1) {
    yd.push_back(update->dt);
  } else {
    cout << "Error: in Plot::compute_dt, x_or_y == " << x_or_y << ". Expected 0 or 1!" << endl;
    exit(1);
  }
}

void Plot::compute_time(string name, int x_or_y)
{
  if (x_or_y == 0) {
    xd.push_back(update->atime);
  } else if (x_or_y == 1) {
    yd.push_back(update->atime);
  } else {
    cout << "Error: in Plot::compute_time, x_or_y == " << x_or_y << ". Expected 0 or 1!" << endl;
    exit(1);
  }
}

void Plot::compute_var(string name, int x_or_y)
{
  if (x_or_y == 0) {
    xd.push_back((*input->vars)[name].result(mpm));
  } else if (x_or_y == 1) {
    yd.push_back((*input->vars)[name].result(mpm));
  } else {
    cout << "Error: in Plot::compute_var, x_or_y == " << x_or_y << ". Expected 0 or 1!" << endl;
    exit(1);
  }
}

void Plot::modify(vector<string> args)
{
  if (args.size() < 1) {
    cout << "Erro: too few arguments given to plot_modify command." << endl;
    exit(1);
  }

  style = args[0];

  vector<string> keyword;

  keyword.push_back(args[2]); // x
  keyword.push_back(args[3]); // y
  clearfields();
  parse_keywords(keyword);
}


void Plot::show()
{
  if ((this->field[0].typeflag == INT) && (this->field[1].typeflag == INT)){
    plt::plot(xi,yi);
  } else if ((this->field[0].typeflag == FLOAT) && (this->field[1].typeflag == INT)){
    plt::plot(xd,yi);
  } else if ((this->field[0].typeflag == BIGINT) && (this->field[1].typeflag == INT)){
    plt::plot(xbi,yi);
  } else if ((this->field[0].typeflag == INT) && (this->field[1].typeflag == FLOAT)){
    plt::plot(xi,yd);
  } else if ((this->field[0].typeflag == FLOAT) && (this->field[1].typeflag == FLOAT)){
    plt::plot(xd,yd);
  } else if ((this->field[0].typeflag == BIGINT) && (this->field[1].typeflag == FLOAT)){
    plt::plot(xbi,yd);
  } else if ((this->field[0].typeflag == INT) && (this->field[1].typeflag == BIGINT)){
    plt::plot(xi,ybi);
  } else if ((this->field[0].typeflag == FLOAT) && (this->field[1].typeflag == BIGINT)){
    plt::plot(xd,ybi);
  } else if ((this->field[0].typeflag == BIGINT) && (this->field[1].typeflag == BIGINT)){
    plt::plot(xbi,ybi);
  }
  plt::xlabel(this->field[0].name);
  plt::ylabel(this->field[1].name);

  string fpng = ("./" + id) + ".png";
  cout << "Save plot as : " << fpng << endl;
  plt::save(fpng);
  plt::show();
}
