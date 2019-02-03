#include "domain.h"
#include "region.h"
#include "style_region.h"

using namespace std;

Domain::Domain(MPM *mpm) : Pointers(mpm)
{
  dimension = 3;

  region_map = new RegionCreatorMap();

#define REGION_CLASS
#define RegionStyle(key,Class) \
  (*region_map)[#key] = &region_creator<Class>;
#include "style_region.h"
#undef RegionStyle
#undef REGION_CLASS
}


Domain::~Domain()
{
  delete region_map;
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

  string *estyle = &args[1];

  if (region_map->find(*estyle) != region_map->end()) {
    RegionCreator region_creator = (*region_map)[*estyle];
    regions.push_back(region_creator(mpm, args));
    regions.back()->init();
    return;
  }
  
}

int Domain::find_region(string name)
{
  for (int iregion = 0; iregion < regions.size(); iregion++)
    if (name.compare(regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   one instance per region style in style_region.h
------------------------------------------------------------------------- */

template <typename T>
Region *Domain::region_creator(MPM *mpm, vector<string> args)
{
  return new T(mpm, args);
}
