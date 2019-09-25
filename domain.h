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
#include <Eigen/Eigen>

using namespace std;

class Domain : protected Pointers {
 public:
  int dimension;                         // 2 = 2d, 3 = 3d

  double boxlo[3],boxhi[3];              // orthogonal box global bounds
  double sublo[3],subhi[3];              // sub-box bounds on this proc

  vector<class Region *> regions;        // list of defined Regions
  vector<class Solid *> solids;          // list of defined Solids

  class Grid *grid;                      // common background grid

  Domain(class MPM *);
  virtual ~Domain();

  void create_domain(vector<string>);
  void set_local_box();
  void add_region(vector<string>);
  int find_region(string);
  void add_solid(vector<string>);
  int find_solid(string);

  typedef Region *(*RegionCreator)(MPM *,vector<string>);
  typedef map<string,RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;

  // typedef Solid *(*SolidCreator)(MPM *,vector<string>);
  // typedef map<string,SolidCreator> SolidCreatorMap;
  // SolidCreatorMap *solid_map;

  int inside(Eigen::Vector3d);

 private:
  template <typename T> static Region *region_creator(MPM *,vector<string>);
  // template <typename T> static Solid *solid_creator(MPM *,vector<string>);
};

#endif
