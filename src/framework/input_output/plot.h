/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#ifndef MPM_PLOT_H
#define MPM_PLOT_H

#include <pointers.h>
#include <vector>

class Plot;

typedef void (Plot::*FnPtrPlt)(string, int);

struct FieldPlt {
  string name = "";
  FnPtrPlt vfunc = nullptr;
  int typeflag = -1;
  int x_or_y = 0;
};

class Plot : protected Pointers {
public:
  string id;
  string style;
  vector<FieldPlt> field; 

  int ivalue;            // integer value to print
  float dvalue;         // float value to print
  bigint bivalue;        // big integer value to print
  vector<float> x;      // vector of values to plot as x
  vector<float> y;      // vector of values to plot as y

  Plot(class MPM*, vector<string>);
  ~Plot() {};

  void modify(vector<string>);

  void parse_keywords(vector<string>);
  void addfield(string, FnPtrPlt, int, int);
  void clearfields();
  void init();
  void save();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step(string, int); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_dt(string, int);   // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_time(string, int); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_var(string, int);
};

#endif
