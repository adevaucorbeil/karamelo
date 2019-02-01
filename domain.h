/* -*- c++ -*- ----------------------------------------------------------*/


#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include <math.h>
#include "pointers.h"
#include "region.h"
#include <map>
#include <string>

using namespace std;

class Domain : protected Pointers {
 public:
  int dimension;                         // 2 = 2d, 3 = 3d

  Domain(class MPM *);
  virtual ~Domain();

  typedef Region *(*RegionCreator)(MPM *,string *);
  typedef map<string,RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;

 private:
  template <typename T> static Region *region_creator(MPM *,string *);
};

#endif
