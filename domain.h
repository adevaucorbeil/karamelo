/* -*- c++ -*- ----------------------------------------------------------*/


#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include <math.h>
#include "pointers.h"
#include "region.h"
#include "solid.h"
#include <map>
#include <string>
#include <vector>

using namespace std;

class Domain : protected Pointers {
 public:
  int dimension;                         // 2 = 2d, 3 = 3d
  vector<class Region *> regions;        // list of defined Regions
  vector<class Solid *> solids;          // list of defined Solids

  Domain(class MPM *);
  virtual ~Domain();

  void add_region(vector<string>);
  int find_region(string);
  void add_solid(vector<string>);
  int find_solid(string);

  typedef Region *(*RegionCreator)(MPM *,vector<string>);
  typedef map<string,RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;

  typedef Solid *(*SolidCreator)(MPM *,vector<string>);
  typedef map<string,SolidCreator> SolidCreatorMap;
  SolidCreatorMap *solid_map;



 private:
  template <typename T> static Region *region_creator(MPM *,vector<string>);
  template <typename T> static Solid *solid_creator(MPM *,vector<string>);
};

#endif
