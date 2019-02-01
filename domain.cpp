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
}


/* ----------------------------------------------------------------------
   one instance per region style in style_region.h
------------------------------------------------------------------------- */

template <typename T>
Region *Domain::region_creator(MPM *mpm, string * args)
{
  return new T(mpm, args);
}
