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

/*! Stores the dimension and geometrical aspects of the simulation box.
 *
 * Stores the simulation dimensions (i.e. 1D, 2D or 3D), the simulation box geometry, 
 * the list of the user defined geometric regions and solids, as well as the background 
 * grid if updated Lagrangian is used.
 */
class Domain : protected Pointers {
 public:
  int dimension;                         ///< 1 == 1D, 2 == 2d, 3 == 3d
  bool created;                          ///< has the domain been created?
  bool axisymmetric;                     ///< true if axisymmetric,  false otherwise
  tagint np_total;                       ///< total number of particles

  double boxlo[3];                       ///< Lower orthogonal box global bounds
  double boxhi[3];                       ///< Higher orthogonal box global bounds
  double sublo[3];                       ///< Lower sub-box bounds on this proc
  double subhi[3];                       ///< Higher sub-box bounds on this proc

  vector<class Region *> regions;        ///< list of defined Regions
  vector<class Solid *> solids;          ///< list of defined Solids

  class Grid *grid;                      ///< pointer to the common background grid (if using Updated Lagrangian)

  Domain(class MPM *);                   ///< Creator.
  virtual ~Domain();                     ///< Destructor.

  void create_domain(vector<string>);    ///< Deprecated function
  void set_dimension(vector<string>);    ///< Called when user calls dimension()
  void set_axisymmetric(vector<string>); ///< Called when user calls axisymmetric()
  bool inside_subdomain(double, double, double); ///< Checks if the set of coordinates lies in the simulation domain.
  bool inside_subdomain_extended(double, double, double, double); ///< Checks if the set of coordinates lies in this proc sub-domain.
  void set_local_box();                  ///< Determine the boundaries of this proc subdomain
  void add_region(vector<string>);       ///< Create a new region
  int find_region(string);               ///< Finds the ID of a region
  void add_solid(vector<string>);        ///< Create a new solid
  int find_solid(string);                ///< Finds the ID of a solid

  typedef Region *(*RegionCreator)(MPM *, vector<string>);
  typedef map<string, RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;          ///< Map of all the known region types

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

/*! \defgroup dimension Dimension
  
\section Usage Usage
dimension(1, domain xmin, domain xmax, cell size)
*/
