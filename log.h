/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_LOG_H
#define MPM_LOG_H

#include "pointers.h"
#include <vector>

class Log;

typedef void (Log::*FnPtr)();

struct Field {
  string name = "";
  FnPtr vfunc = NULL;
  int typeflag;
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

  void parse_keywords(vector<string>);
  void addfield(string, FnPtr, int);
  void header();
  void init();

  void write();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();
  void compute_dt();
  void compute_time();
};

#endif
