/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_PLOT_H
#define MPM_PLOT_H

#include "pointers.h"
#include <vector>

class Plot;

typedef void (Plot::*FnPtrPlt)(string, int);

struct FieldPlt {
  string name = "";
  FnPtrPlt vfunc = NULL;
  int typeflag = -1;
  int x_or_y = 0;
};

class Plot : protected Pointers {
public:
  string id;
  string style;
  vector<FieldPlt> field; 

  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print
  vector<int> xi;        // vector of integer values to plot as x
  vector<int> yi;        // vector of integer values to plot as y
  vector<double> xd;     // vector of double values to plot as x
  vector<double> yd;     // vector of double values to plot as y
  vector<bigint> xbi;    // vector of big integer values to plot as x
  vector<bigint> ybi;    // vector of big integer values to plot as y

  Plot(class MPM*, vector<string>);
  ~Plot() {};

  void modify(vector<string>);

  void parse_keywords(vector<string>);
  void addfield(string, FnPtrPlt, int, int);
  void clearfields();
  void init();
  void save();
  void show();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step(string, int); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_dt(string, int);   // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_time(string, int); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_var(string, int);
};

#endif
