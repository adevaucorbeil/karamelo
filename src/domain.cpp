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


/*! Initializes variables and generates the map of known region types (at compile time).
 */
Domain::Domain(MPM *mpm) : Pointers(mpm)
{
  dimension    = 3;
  created      = false;
  axisymmetric = false;
  np_total     = 0;
  np_local     = 0;

  boxlo[0] = boxlo[1] = boxlo[2] = 0;
  boxhi[0] = boxhi[1] = boxhi[2] = 0;

  region_map = new RegionCreatorMap();
  // solid_map = new SolidCreatorMap();


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

/*! Destroys the lists of regions and solids (Domain::regions, Domain::solids) 
 * as well as the map of known region types.
 */
Domain::~Domain()
{
  for (int i = 0; i < regions.size(); i++)
    delete regions[i];
  for (int i = 0; i < solids.size(); i++)
    delete solids[i];

  delete region_map;
  // delete solid_map;

  if (!update->method->is_TL) delete grid;
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

/*! This function is the C++ equivalent to the region() user function.\n
 * Syntax: region(name, type, type specific arguments)\n
 * This function checks if the region name was already used.\n
 * If not, it creates an entry in the vector Domain::regions and calls Region::Region()
 */
void Domain::add_region(vector<string> args)
{
  // cout << "In add_region" << endl;

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

/*! This function checks if 'name' is already used for a region.\n
 * If a region named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
int Domain::find_region(string name)
{
  if (regions.size() > 0)
  {
    for (int iregion = 0; iregion < regions.size(); iregion++)
    {
      //cout << "regions[" << iregion << "]->id=" << regions[iregion]->id << endl;
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

/*! This function is the C++ equivalent to the solid() user function.\n
 * This function checks if the solid name was already used.\n
 * If not, it creates an entry in the vector Domain::solids.
 */
void Domain::add_solid(vector<string> args)
{
  // cout << "In add_solid" << endl;

  if (find_solid(args[0]) >= 0) error->all(FLERR, "Error: reuse of solid ID.\n");

  // create the Solid
  solids.push_back(new Solid(mpm, args));
  solids.back()->init();
}

/*! This function checks if 'name' is already used for a solid.\n
 * If a solid named 'name' exists, it returns its ID. It returns -1 otherwise.
 */
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

/*! inside = 1 if x,y,z is inside or on the boundary of the simulation domain.
 *  inside = 0 if x,y,z is outside and not on boundary of the simulation domain.
 */
int Domain::inside(Eigen::Vector3d x)
{
  // cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside
  // the domain" << endl;
  if (x[0] >= boxlo[0] && x[0] <= boxhi[0] && x[1] >= boxlo[1] &&
      x[1] <= boxhi[1] && x[2] >= boxlo[2] && x[2] <= boxhi[2])
    return 1;
  return 0;
}

/*! inside = 1 if x,y,z is inside or on the boundary of this proc domain.
 *  inside = 0 if x,y,z is outside and not on boundary of this proc domain.
 */
bool Domain::inside_subdomain(double x, double y, double z) {
  if (x < sublo[0]) return false;
  if (x > subhi[0]) return false;
  if (y < sublo[1]) return false;
  if (y > subhi[1]) return false;
  if (z < sublo[2]) return false;
  if (z > subhi[2]) return false;
  return true;
}

/*! Determine the CPU owning the particle with x, y, z coordinates.
 */
int Domain::which_CPU_owns_me(double x, double y, double z) {
  if (x < sublo[0]) {
    // target[0] = -1;
    if (y < sublo[1]) {
      // target[1] = -1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {-1, -1, -1};
	return universe->procneigh[2][0] - 1 - universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {-1, -1, 1};
	return universe->procneigh[2][1] - 1 - universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {-1, -1, 0};
	return universe->procneigh[1][0] - 1;
      }
    } else if (y > subhi[1]) {
      // target[1] = 1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {-1, 1, -1};
	return universe->procneigh[2][1] - 1 + universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {-1, 1, 1};
	return universe->procneigh[2][1] - 1 + universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {-1, 1, 0};
	return universe->procneigh[1][1] - 1;
      }
    } else {
      // target[1] = 0;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {-1, 0, -1};
	return universe->procneigh[2][0] - 1;
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {-1, 0, 1};
	return universe->procneigh[2][1] - 1;
      } else {
        // target[2] = 0;
	// target = {-1, 0, 0};
	return universe->procneigh[0][0];
      }
    }
  } else if (x > subhi[0]) {
    // target[0] = 1;
    if (y < sublo[1]) {
      // target[1] = -1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {1, -1, -1};
	return universe->procneigh[2][0] + 1 - universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {1, -1, 1};
	return universe->procneigh[2][1] + 1 - universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {1, -1, 0};
	return universe->procneigh[1][0] + 1;
      }
    } else if (y > subhi[1]) {
      // target[1] = 1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {1, 1, -1};
	return universe->procneigh[2][0] + 1 + universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {1, 1, 1};
	return universe->procneigh[2][1] + 1 + universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {1, 1, 0};
	return universe->procneigh[1][1] + 1;
      }
    } else {
      // target[1] = 0;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {1, 0, -1};
	return universe->procneigh[2][0] + 1;
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {1, 0, 1};
	return universe->procneigh[2][1] + 1;
      } else {
        // target[2] = 0;
	// target = {1, 0, 0};
	return universe->procneigh[0][1];
      }
    }
  } else {
    // target[0] = 0;
    if (y < sublo[1]) {
      // target[1] = -1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {0, -1, -1};
	return universe->procneigh[2][0] - universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {0, -1, 1};
	return universe->procneigh[2][1] - universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {0, -1, 0};
	return universe->procneigh[1][0];
      }
    } else if (y > subhi[1]) {
      // target[1] = 1;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {0, 1, -1};
	return universe->procneigh[2][0] + universe->procgrid[1];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {0, 1, 1};
	return universe->procneigh[2][1] + universe->procgrid[1];
      } else {
        // target[2] = 0;
	// target = {0, 1, 0};
	return universe->procneigh[1][1];
      }
    } else {
      // target[1] = 0;
      if (z < sublo[2]) {
        // target[2] = -1;
	// target = {0, 0, -1};
	return universe->procneigh[2][0];
      } else if (z > subhi[2]) {
        // target[2] = 1;
	// target = {0, 0, 1};
	return universe->procneigh[2][1];
      } else {
        // target[2] = 0;
	// target = {0, 0, 0};
	return universe->me;
      }
    }
  }

  return -1;
}


/*! inside = 1 if x,y,z is inside or on the boundary of this proc domain, extended by h.
 *  inside = 0 if x,y,z is outside and not on boundary of this proc domain, extended by h.
 */
bool Domain::inside_subdomain_extended(double x, double y, double z, double h) {
  if (x < sublo[0] - h) return false;
  if (x > subhi[0] + h) return false;
  if (y < sublo[1] - h) return false;
  if (y > subhi[1] + h) return false;
  if (z < sublo[2] - h) return false;
  if (z > subhi[2] + h) return false;
  return true;
}

/*! The boundaries of the domain attributed to this proc are stored in 
 * Domain::sublo and Domain::subhi variables.
 */
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

/*! This function is the C++ equivalent to the dimension() user function.\n
 * Syntax: dimension(N, xmin, xmax, ymin, ymax, zmin, zmax, h)
 * N: 1 for 1D, 2 for 2D, 3 for 3D\.n
 * xmin, xmax: domain boundaries in x.\n
 * ymin, ymax: domain boundaries in y -- Optional when N == 1.\n
 * zmin, zmax: domain boundaries in z -- Optional when N == 2.\n
 * h: background grid cell size. Not optional if using Total Lagrangian, but value not used.
 */
void Domain::set_dimension(vector<string> args) {

  if (args.size() == 0) {
    string error_str = "Error: dimension did not receive enough arguments\n";
    for (auto &x : usage_dimension)
      error_str += x.second;
    error->all(FLERR, error_str);
  }

  if (!update->method->is_TL) grid = new Grid(mpm);

  int m = 0;
  int dim = (int)input->parsev(args[m]);

  if (dim != 1 && dim != 2 && dim != 3) {
    error->all(FLERR, "Error: dimension argument: " + args[m] + "\n.");
  } else {
    dimension = dim;
  }

  if (universe->me == 0)
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


/*! This function is the C++ equivalent to the axisymmetric() user function.\n
 * Syntax: axisymmetric(true) or axisymmetric(false)\n
 * Sets Domain::axisymmetric to true or false. Note that this variable is initialized as false.
 */
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

/*! Write box bounds, list of regions and solids to restart file
 */
void Domain::write_restart(ofstream* of){

  // Write boxlo:
  of->write(reinterpret_cast<const char *>(&boxlo[0]), 3*sizeof(double));
  
  // Write boxhi:
  of->write(reinterpret_cast<const char *>(&boxhi[0]), 3*sizeof(double));

  // Write sublo:
  of->write(reinterpret_cast<const char *>(&sublo[0]), 3*sizeof(double));

  // Write subhi:
  of->write(reinterpret_cast<const char *>(&subhi[0]), 3*sizeof(double));

  // Write axisymmetric:
  of->write(reinterpret_cast<const char *>(&axisymmetric), sizeof(bool));

  // Write  np_total:
  of->write(reinterpret_cast<const char *>(&np_total), sizeof(tagint));
  // cout << "np_total=" << np_total << endl;
  
  if (!update->method->is_TL) {
    of->write(reinterpret_cast<const char *>(&grid->cellsize), sizeof(double));
    // cout << "cellsize=" << grid->cellsize << endl;
  }

  // Save regions:
  size_t N = regions.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Nr = regions[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Nr), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(regions[i]->id.c_str()), Nr);
    // cout << "id = " << regions[i]->id << endl;

    Nr = regions[i]->style.size();
    of->write(reinterpret_cast<const char *>(&Nr), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(regions[i]->style.c_str()), Nr);
    regions[i]->write_restart(of);
    // cout << "style = " << regions[i]->style << endl;
  }

  // Save materials:
  material->write_restart(of);

  // Save solids:
  N = solids.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(int));

  for (int i = 0; i < N; i++) {
    size_t Ns = solids[i]->id.size();
    of->write(reinterpret_cast<const char *>(&Ns), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(solids[i]->id.c_str()), Ns);
    // cout << "id = " << solids[i]->id << endl;
    solids[i]->write_restart(of);
  }
}

/*! Read box bounds, list of regions and solids to restart file
 */
void Domain::read_restart(ifstream *ifr) {
  // Write boxlo:
  ifr->read(reinterpret_cast<char *>(&boxlo[0]), 3*sizeof(double));
  // cout << "boxlo=[" << boxlo[0] << "," << boxlo[1] << "," << boxlo[2] << endl;
  
  // Write boxhi:
  ifr->read(reinterpret_cast<char *>(&boxhi[0]), 3*sizeof(double));
  // cout << "boxhi=[" << boxhi[0] << "," << boxhi[1] << "," << boxhi[2] << endl;

  // Write sublo:
  ifr->read(reinterpret_cast<char *>(&sublo[0]), 3*sizeof(double));
  // cout << "sublo=[" << sublo[0] << "," << sublo[1] << "," << sublo[2] << endl;

  // Write subhi:
  ifr->read(reinterpret_cast<char *>(&subhi[0]), 3*sizeof(double));
  // cout << "subhi=[" << subhi[0] << "," << subhi[1] << "," << subhi[2] << endl;

  // Write axisymmetric:
  ifr->read(reinterpret_cast<char *>(&axisymmetric), sizeof(bool));
  // cout << "axisymmetric=" << axisymmetric << endl;

  // Write  np_total:
  ifr->read(reinterpret_cast<char *>(&np_total), sizeof(tagint));
  // cout << "np_total=" << np_total << endl;
  
  universe->set_proc_grid();
  if (!update->method->is_TL) {
    grid = new Grid(mpm);
    ifr->read(reinterpret_cast<char *>(&grid->cellsize), sizeof(double));
    // cout << "cellsize=" << grid->cellsize << endl;
    
    grid->init(boxlo, boxhi);
  }
  created = true;


  // Pull regions:
  size_t N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  regions.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Nr = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Nr), sizeof(size_t));
    id.resize(Nr);

    ifr->read(reinterpret_cast<char *>(&id[0]), Nr);
    // cout << "id = " << id << endl;

    string style = "";
    ifr->read(reinterpret_cast<char *>(&Nr), sizeof(size_t));
    style.resize(Nr);

    ifr->read(reinterpret_cast<char *>(&style[0]), Nr);
    // cout << "style = " << style << endl;
    RegionCreator region_creator = (*region_map)[style];
    regions[i] = region_creator(mpm, vector<string>{id, style, "restart"});
    regions[i]->read_restart(ifr);
    regions[i]->init();
  }

  // Save materials:
  material->read_restart(ifr);

  // Read solids:
  N = 0;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(int));
  solids.resize(N);

  for (int i = 0; i < N; i++) {
    size_t Ns = 0;
    string id = "";

    ifr->read(reinterpret_cast<char *>(&Ns), sizeof(size_t));
    id.resize(Ns);

    ifr->read(reinterpret_cast<char *>(&id[0]), Ns);
    // cout << "id = " << id << endl;

    solids[i] = new Solid(mpm, vector<string>{id, "restart"});
    solids[i]->read_restart(ifr);
    solids[i]->init();
  }
}
