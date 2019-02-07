/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_OUTPUT_H
#define MPM_OUTPUT_H

#include "pointers.h"
#include <vector>

class Output : protected Pointers {
public:
  bigint next;                 // next timestep for any kind of output

  
  bigint next_log;             // next timestep for log output
  int every_log;               // freq for log output
  class Log *log;

  bigint next_dump_any;        // next timestep for any Dump
  vector<int> every_dump;      // write freq for each Dump, 0 if var
  vector<bigint> next_dump;    // next timestep to do each Dump
  vector<bigint> last_dump;    // last timestep each snapshot was output
  int ndumps;                  // number of defined Dumps, should always be equal to dumps.size()
  vector<class Dump *> dumps;  // list of defined Dumps

  Output(class MPM *);
  ~Output();

  void setup();
  void write(bigint);                  // output for current timestep

  void set_log(vector<string>);        // set log output freqquency

  void add_dump(vector<string>);       // add a Dump to Dump list
  int find_dump(string);               // find a Dump in Dump list
  void modify_dump(vector<string>);    // modify a Dump
  void delete_dump(string);            // delete a Dump from Dump list
};

#endif
