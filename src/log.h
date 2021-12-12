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

#ifndef MPM_LOG_H
#define MPM_LOG_H

#include <pointers.h>
#include <vector>

class Log;

typedef void (Log::*FnPtr)(string);

struct Field {
  string name = "";
  FnPtr vfunc = nullptr;
  int typeflag = -1;
};

class Log : protected Pointers {
public:
  string style;
  vector<Field> field; 

  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print

  Log(class MPM*, vector<string>);
  ~Log() {};

  void modify(vector<string>);

  void parse_keywords(vector<string>);
  void addfield(string, FnPtr, int);
  void clearfields();
  void header();
  void init();

  void write();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step(string); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_dt(string);   // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_time(string); // Here the argument is not required, but prototype should be the same for all compute functions.
  void compute_var(string);
};

#endif
