#include "domain.h"
#include "region.h"
#include "memory.h"
#include "style_region.h"
#include <Eigen/Eigen>

using namespace std;

Domain::Domain(MPM *mpm) : Pointers(mpm)
{
  dimension = 3;

  boxlo[0] = boxlo[1] = boxlo[2] = 0;
  boxhi[0] = boxhi[1] = boxhi[2] = 0;

  region_map = new RegionCreatorMap();
  // solid_map = new SolidCreatorMap();

  grid = new Grid(mpm);

#define REGION_CLASS
#define RegionStyle(key,Class) \
  (*region_map)[#key] = &region_creator<Class>;
#include "style_region.h"
#undef RegionStyle
#undef REGION_CLASS

// #define SOLID_CLASS
// #define SolidStyle(key,Class) \
//   (*solid_map)[#key] = &solid_creator<Class>;
// #include "style_solid.h"
// #undef SolidStyle
// #undef SOLID_CLASS
}

Domain::~Domain()
{
  for (int i = 0; i < regions.size(); i++) delete regions[i];
  for (int i = 0; i < solids.size(); i++) delete solids[i];

  delete region_map;
  // delete solid_map;

  delete grid;
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Domain::add_region(vector<string> args){
  cout << "In add_region" << endl;

  if (find_region(args[0]) >= 0) {
    cout << "Error: reuse of region ID" << endl;
    exit(1);
  }

    // create the Region

  if (region_map->find(args[1]) != region_map->end()) {
    RegionCreator region_creator = (*region_map)[args[1]];
    regions.push_back(region_creator(mpm, args));
    regions.back()->init();

  }
  else {
    cout << "Unknown region style " << args[1] << endl;
    exit(1);
  }
  

}

int Domain::find_region(string name)
{
  cout << "regions.size = " << regions.size() << endl;
  if (regions.size()>0) {
    for (int iregion = 0; iregion < regions.size(); iregion++) {
      cout << "regions["<< iregion <<"]->id=" << regions[iregion]->id << endl;
      if (name.compare(regions[iregion]->id) == 0) return iregion;
    }
    return -1;
  } else return -1;
}

/* ----------------------------------------------------------------------
   one instance per region style in style_region.h
------------------------------------------------------------------------- */

template <typename T>
Region *Domain::region_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}

/* ----------------------------------------------------------------------
   create a new solid
------------------------------------------------------------------------- */

void Domain::add_solid(vector<string> args){
  cout << "In add_solid" << endl;

  if (find_solid(args[0]) >= 0) {
    cout << "Error: reuse of solid ID" << endl;
    exit(1);
  }

  // create the Solid

  // string *estyle = &args[1];

  // if (solid_map->find(*estyle) != solid_map->end()) {
  //   SolidCreator solid_creator = (*solid_map)[*estyle];
  //   solids.push_back(solid_creator(mpm, args));
  //   solids.back()->init();
  // }
  // else {
  //   cout << "Unknown solid style " << *estyle << endl;
  //   exit(1);
  // }
  solids.push_back(new Solid(mpm, args));
  solids.back()->init();
  
}

int Domain::find_solid(string name)
{
  for (int isolid = 0; isolid < solids.size(); isolid++)
    if (name.compare(solids[isolid]->id) == 0) return isolid;
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per solid style in style_solid.h
------------------------------------------------------------------------- */

// template <typename T>
// Solid *Domain::solid_creator(MPM *mpm, vector<string> args)
// {
//   return new T(mpm, args);
// }


/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Domain::inside(Eigen::Vector3d x)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the domain" << endl;
  if (   x[0] >= boxlo[0] && x[0] <= boxhi[0]
      && x[1] >= boxlo[1] && x[1] <= boxhi[1]
      && x[2] >= boxlo[2] && x[2] <= boxhi[2])
    return 1;
  return 0;
}

void Domain::set_local_box() {
  if (dimension == 1) {
    ;
  } else if (dimension == 2) {
    ;
  } else {
    ;
  }
}
