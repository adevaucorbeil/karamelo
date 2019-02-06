/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_OUTPUT_H
#define MPM_OUTPUT_H

#include "pointers.h"
#include <vector>

class Output : protected Pointers {
public:
  bigint next_dump_any;        // next timestep for any Dump
  vector<int> every_dump;      // write freq for each Dump, 0 if var
  bigint *next_dump;           // next timestep to do each Dump
  bigint *last_dump;           // last timestep each snapshot was output
  char **var_dump;             // variable name for dump frequency
  int *ivar_dump;              // variable index for dump frequency
  vector<class Dump *> dumps;  // list of defined Dumps

  Output(class MPM *);
  ~Output();

  void add_dump(vector<string>);       // add a Dump to Dump list
  int find_dump(string);               // find a Dump in Dump list
  void modify_dump(vector<string>);    // modify a Dump
  void delete_dump(string);            // delete a Dump from Dump list
};

#endif
