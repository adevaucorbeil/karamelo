/* ----------------------------------------------------------------------
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

#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "method.h"
#include "region.h"
#include "style_region.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>

using namespace std;

#ifdef DEBUG
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
#endif

Domain::Domain(MPM *mpm) : Pointers(mpm)
{
  dimension    = 3;
  created      = false;
  axisymmetric = false;
  np_total     = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = 0;
  boxhi[0] = boxhi[1] = boxhi[2] = 0;

  region_map = new RegionCreatorMap();
  // solid_map = new SolidCreatorMap();

#ifdef DEBUG
  plt::figure_size(1200, 780);
#endif

#define REGION_CLASS
#define RegionStyle(key, Class) (*region_map)[#key] = &region_creator<Class>;
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
  for (int i = 0; i < regions.size(); i++)
    delete regions[i];
  for (int i = 0; i < solids.size(); i++)
    delete solids[i];

  delete region_map;
  // delete solid_map;

  delete grid;
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Domain::add_region(vector<string> args)
{
  cout << "In add_region" << endl;

  if (find_region(args[0]) >= 0) error->all(FLERR, "Error: reuse of region ID.\n");

  // create the Region

  if (region_map->find(args[1]) != region_map->end())
  {
    RegionCreator region_creator = (*region_map)[args[1]];
    regions.push_back(region_creator(mpm, args));
    regions.back()->init();
  }
  else error->all(FLERR, "Unknown region style " + args[1] + "\n");
}

int Domain::find_region(string name)
{
  cout << "regions.size = " << regions.size() << endl;
  if (regions.size() > 0)
  {
    for (int iregion = 0; iregion < regions.size(); iregion++)
    {
      cout << "regions[" << iregion << "]->id=" << regions[iregion]->id << endl;
      if (name.compare(regions[iregion]->id) == 0)
        return iregion;
    }
    return -1;
  }
  else
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

/* ----------------------------------------------------------------------
   create a new solid
------------------------------------------------------------------------- */

void Domain::add_solid(vector<string> args)
{
  cout << "In add_solid" << endl;

  if (find_solid(args[0]) >= 0) error->all(FLERR, "Error: reuse of region ID.\n");

  // create the Solid
  solids.push_back(new Solid(mpm, args));
  solids.back()->init();
}

int Domain::find_solid(string name)
{
  for (int isolid = 0; isolid < solids.size(); isolid++)
    if (name.compare(solids[isolid]->id) == 0)
      return isolid;
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
  // cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside
  // the domain" << endl;
  if (x[0] >= boxlo[0] && x[0] <= boxhi[0] && x[1] >= boxlo[1] &&
      x[1] <= boxhi[1] && x[2] >= boxlo[2] && x[2] <= boxhi[2])
    return 1;
  return 0;
}

bool Domain::inside_subdomain(double x, double y, double z) {
  if (x < sublo[0]) return false;
  if (x > subhi[0]) return false;
  if (y < sublo[1]) return false;
  if (y > subhi[1]) return false;
  if (z < sublo[2]) return false;
  if (z > subhi[2]) return false;
  return true;
}

bool Domain::inside_subdomain_extended(double x, double y, double z, double h) {
  if (x < sublo[0] - h) return false;
  if (x > subhi[0] + h) return false;
  if (y < sublo[1] - h) return false;
  if (y > subhi[1] + h) return false;
  if (z < sublo[2] - h) return false;
  if (z > subhi[2] + h) return false;
  return true;
}

void Domain::set_local_box() {
  int *procgrid = universe->procgrid;
  int *myloc = universe->myloc;

  double l[3] = {boxhi[0] - boxlo[0],
		 boxhi[1] - boxlo[1],
		 boxhi[2] - boxlo[2]};

  double h[3];
  h[0] = l[0]/procgrid[0];
  sublo[0] = myloc[0]*h[0] + boxlo[0];
  subhi[0] = sublo[0] + h[0];

  if (dimension >= 2) {
    h[1] = l[1]/procgrid[1];
    sublo[1] = myloc[1]*h[1] + boxlo[1];
    subhi[1] = sublo[1] + h[1];
  } else {
    sublo[1] = subhi[1] = 0;
  }

  if (dimension == 3) {
    h[2] = l[2]/procgrid[2];
    sublo[2] = myloc[2]*h[2] + boxlo[2];
    subhi[2] = sublo[2] + h[2];
  } else {
    sublo[2] = subhi[2] = 0;
  }

#ifdef DEBUG
  cout << "proc " << universe->me << "\tboxlo=[" << boxlo[0] << "," << boxlo[1] << "," << boxlo[2] << "]\t boxhi=["<< boxhi[0] << "," << boxhi[1] << "," << boxhi[2] << "]\n";
  cout << "proc " << universe->me << "\tsublo=[" << sublo[0] << "," << sublo[1] << "," << sublo[2] << "]\t subhi=["<< subhi[0] << "," << subhi[1] << "," << subhi[2] << "]\n";
#endif
}

void Domain::create_domain(vector<string> args) {

  // Check that a method is available:
  if (update->method == NULL)
    error->all(FLERR, "Error: a method should be defined before calling create_domain()!\n");

  if (!update->method->is_TL) grid = new Grid(mpm);

  if (dimension ==1) {
    if (update->method->is_TL) {
      if (args.size() < 2)
	error->all(FLERR, "Error: create_domain received not enough arguments: 2 needed for 1D total Lagrangian simulations: domain xmin, domain xmax.\n");
      else if (args.size() > 2)
	error->all(FLERR, "Error: create_domain received too many arguments: 2 needed for 1D total Lagrangian simulations: domain xmin, domain xmax.\n");
    } else {
      if (args.size() < 3)
	error->all(FLERR, "Error: create_domain received not enough arguments: 3 needed for 1D total Lagrangian simulations: domain xmin, domain xmax, grid cell size.\n");
      else if (args.size() > 3)
	error->all(FLERR, "Error: create_domain received too many arguments: 3 needed for 1D total Lagrangian simulations: domain xmin, domain xmax, grid cell size.\n");
    }

    boxlo[0] = (double) input->parsev(args[0]);
    boxhi[0] = (double) input->parsev(args[1]);

    if (!update->method->is_TL) {
      grid->cellsize = (double) input->parsev(args[2]);

      if (grid->cellsize < 0) 
	error->all(FLERR, "Error: cellsize negative! You gave: " + to_string(grid->cellsize) + "\n");

      grid->init(boxlo, boxhi);
    }
  } else if (dimension == 2) {
    if (update->method->is_TL) {
      if (args.size() < 4)
	error->all(FLERR, "Error: create_domain received not enough arguments: 2 needed for 2D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax.\n");
      else if (args.size() > 4)
	error->all(FLERR, "Error: create_domain received too many arguments: 2 needed for 2D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax.\n");
    } else {
      if (args.size() < 5)  
	error->all(FLERR, "Error: create_domain received not enough arguments: 3 needed for 2D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, grid cell size.\n");
      else if (args.size() > 5)  
	error->all(FLERR, "Error: create_domain received too many arguments: 3 needed for 2D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, grid cell size.\n");
    }

    boxlo[0] = (double) input->parsev(args[0]);
    boxhi[0] = (double) input->parsev(args[1]);
    boxlo[1] = (double) input->parsev(args[2]);
    boxhi[1] = (double) input->parsev(args[3]);

    if (!update->method->is_TL) {
      grid->cellsize = (double) input->parsev(args[4]);

      if (grid->cellsize < 0) 
	error->all(FLERR, "Error: cellsize negative! You gave: " + to_string(grid->cellsize) + "\n");

      grid->init(boxlo, boxhi);
    }
  } else {// dimension ==3
    if (update->method->is_TL) {
      if (args.size() < 6)
	error->all(FLERR, "Error: create_domain received not enough arguments: 2 needed for 3D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, domain zmin, domain zmax.\n");
      else if (args.size() > 6)
	error->all(FLERR, "Error: create_domain received too many arguments: 2 needed for 3D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, domain zmin, domain zmax.\n");
    } else {
      if (args.size() < 7)
	error->all(FLERR, "Error: create_domain received not enough arguments: 3 needed for 3D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, domain zmin, domain zmax, grid cell size.\n");
      else if (args.size() > 7)
	error->all(FLERR, "Error: create_domain received too many arguments: 3 needed for 3D total Lagrangian simulations: domain xmin, domain xmax, domain ymin, domain ymax, domain zmin, domain zmax, grid cell size.\n");
    }
    boxlo[0] = (double) input->parsev(args[0]);
    boxhi[0] = (double) input->parsev(args[1]);
    boxlo[1] = (double) input->parsev(args[2]);
    boxhi[1] = (double) input->parsev(args[3]);
    boxlo[2] = (double) input->parsev(args[4]);
    boxhi[2] = (double) input->parsev(args[5]);

    if (!update->method->is_TL) {
      grid->cellsize = (double) input->parsev(args[6]);

      if (grid->cellsize < 0) 
	error->all(FLERR, "Error: cellsize negative! You gave: " + to_string(grid->cellsize) + "\n");

      grid->init(boxlo, boxhi);
    }
  }
  
  // Set proc grid
  universe->set_proc_grid();
  created = true;
}

void Domain::set_dimension(vector<string> args) {

  if (args.size() == 0) {
    string error_str = "Error: dimension did not receive enough arguments\n";
    for (auto &x : usage_dimension)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  int m = 0;
  int dim = (int)input->parsev(args[m]);

  if (dim != 1 && dim != 2 && dim != 3) {
    error->all(FLERR, "Error: dimension argument: " + args[m] + "\n.");
  } else {
    dimension = dim;
  }

  cout << "Set dimension to " << dim << endl;

  if (args.size() < Nargs_dimension.find(args[m])->second) {
    error->all(FLERR, "Error: not enough arguments.\n"
	       + usage_dimension.find(args[m])->second);
  } else if (args.size() > Nargs_dimension.find(args[m])->second) {
    error->all(FLERR, "Error: too many arguments.\n"
	       + usage_dimension.find(args[m])->second);
  }

  boxlo[0] = (double)input->parsev(args[++m]);
  boxhi[0] = (double)input->parsev(args[++m]);
  if (dim > 1) {
    boxlo[1] = (double)input->parsev(args[++m]);
    boxhi[1] = (double)input->parsev(args[++m]);
  }
  if (dim == 3) {
    boxlo[2] = (double)input->parsev(args[++m]);
    boxhi[2] = (double)input->parsev(args[++m]);
  }

  if (!update->method->is_TL) {
    grid->cellsize = (double)input->parsev(args[++m]);

    if (grid->cellsize < 0) {
      error->all(FLERR, "Error: cellsize negative! You gave: " + to_string(grid->cellsize) + "\n.");
    }
  }

  universe->set_proc_grid();
  if (!update->method->is_TL) {
    grid->init(boxlo, boxhi);
  }
  created = true;
}

void Domain::set_axisymmetric(vector<string> args) {
  if (args.size() < 1) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage_axisymmetric);
  } else if (args.size() > 1) {
    error->all(FLERR, "Error: too many arguments.\n" + usage_axisymmetric);
  }

  if (args[0].compare("true") == 0) {
    axisymmetric = true;
  } else if (args[0].compare("false") == 0) {
    axisymmetric = false;
  }
}
