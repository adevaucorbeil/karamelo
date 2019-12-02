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

#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H

#include "pointers.h"
#include "region.h"
#include "solid.h"
#include <Eigen/Eigen>
#include <map>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

class Domain : protected Pointers {
 public:
  int dimension;                         // 2 = 2d, 3 = 3d
  bool created;                          // has the domain been created?
  bool axisymmetric;                     // true or false

  double boxlo[3],boxhi[3];              // orthogonal box global bounds
  double sublo[3],subhi[3];              // sub-box bounds on this proc

  vector<class Region *> regions; // list of defined Regions
  vector<class Solid *> solids;   // list of defined Solids

  class Grid *grid; // common background grid

  Domain(class MPM *);
  virtual ~Domain();

  void create_domain(vector<string>);
  void set_dimension(vector<string>);
  void set_axisymmetric(vector<string>);
  bool inside_subdomain(double, double, double);
  bool inside_subdomain_extended(double, double, double, double);
  void set_local_box();
  void add_region(vector<string>);
  int find_region(string);
  void add_solid(vector<string>);
  int find_solid(string);

  typedef Region *(*RegionCreator)(MPM *, vector<string>);
  typedef map<string, RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;

  // typedef Solid *(*SolidCreator)(MPM *,vector<string>);
  // typedef map<string,SolidCreator> SolidCreatorMap;
  // SolidCreatorMap *solid_map;

  int inside(Eigen::Vector3d);

private:
  template <typename T> static Region *region_creator(MPM *, vector<string>);
  // template <typename T> static Solid *solid_creator(MPM *,vector<string>);

  const map<string, string> usage_dimension = {
      {"1", "Usage: dimension(1, domain xmin, domain xmax, cell size)\n"},
      {"2", "Usage: dimension(2, domain xmin, domain xmax, domain ymin, domain "
            "ymax, cell size)\n"},
      {"3", "Usage: dimension(3, domain xmin, domain xmax, domain ymin, domain "
            "ymax, domain zmin, domain zmax, cell size)\n"}};

  const map<string, int> Nargs_dimension = {
      {"1", 4}, {"2", 6}, {"2_axis", 7}, {"3", 8}};

  string usage_axisymmetric = "axisymmetric(true) or axisymmetric(false)\n";
};

#endif
