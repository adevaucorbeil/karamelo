/* -*- c++ -*- ----------------------------------------------------------*/


#ifndef LMP_GROUP_H
#define LMP_GROUP_H

#include "pointers.h"
#include <string>
#include <vector>

using namespace std;

class Group : protected Pointers {
 public:
  int ngroup;                  // # of defined groups
  string *names;               // name of each group
  int *bitmask;                // one-bit mask for each group
  int *inversemask;            // inverse mask for each group
  string *pon;                 // group of particles if pon == "particles", or nodes if pon = "nodes"
  int *solid;                  // solids corresponding to each group, if -1, all solids are considered.

  Group(class MPM *);
  virtual ~Group();

  void assign(vector<string>); // assign atoms to a new or existing group
  int find(string);            // return group index
  int find_unused();           // return index of first available group

  double xcm(int, int);        // determine the centre of mass of a group
};

#endif
