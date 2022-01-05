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

#include <region_stl.h>
#include <error.h>
#include <domain.h>
#include <universe.h>
#include <update.h>
#include <climits>
#include <cstring>

#include <sstream>

Stl::Stl(MPM *mpm, vector<string> args) : Region(mpm, args) {
  if (universe->me == 0)
    cout << "Initiate Stl" << endl;

  if (args.size() != 3) {
    error->all(FLERR, "Error: wrong number of arguments.\n");
  }

  if (args[2].compare("restart") == 0) {
    // If the keyword restart, we are expecting to have read_restart()
    // launched right after.
    input_file_name = "";
    name = "";
    return;
  }

  input_file_name = args[2];

  cout << "FILE: " << input_file_name.c_str() << endl;

  fstream input_file(input_file_name.c_str(), ios::in | ios::binary);

  if (!input_file) {
    error->all(FLERR, "Error: input file could not be opened.\n");
  }

  char buffer[6];
  input_file.read(buffer, 5);
  buffer[5] = '\0';

  if (!strcmp(buffer, "solid")) { // ascii stl
    string token;
    
    auto require_token = [&](const string &required_token0,
                            const string &required_token1 = "") {                  
      input_file >> token;

      if (token != required_token0 && token != required_token1) {
        stringstream ss;
        ss << "Error: stl file not valid (expected " << required_token0;

        if (!required_token1.empty()) {
          ss << " or " << required_token1;
        }

        ss << " but found " << token << " instead)." << endl;

        error->all(FLERR, ss.str());
      }
    };

    // note: name might be optional?
    input_file >> name;

    while (true) {
      require_token("endsolid", "facet");

      if (token == "endsolid") {
        break;
      }

      facets.emplace_back();
      Facet &facet = facets.back();

      require_token("normal");

      for (int i = 0; i < 3; i++) {
        input_file >> facet.normal(i);
      }

      require_token("outer");
      require_token("loop");

      for (int i = 0; i < 3; i++) {
        require_token("vertex");

        for (int j = 0; j < 3; j++) {
          input_file >> facet.at(i)(j);
        }

        add(facet.at(i));
      }

      require_token("endloop");
      require_token("endfacet");
    }

    require_token(name);
  }
  else { // binary stl
    input_file.ignore(75);
    
    uint32_t n;
    input_file.read(reinterpret_cast<char *>(&n), sizeof n);

    if (CHAR_BIT*sizeof(float) != 32) {
      error->all(FLERR, "Error: floating type not IEEE compliant.\n");
    }

    float value;

    for (int i = 0; i < n; i++) {
      facets.emplace_back();
      Facet &facet = facets.back();

      for (int j = 0; j < 3; j++) {
        input_file.read(reinterpret_cast<char *>(&value), sizeof(float));
        facet.normal(j) = value;
      }

      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          input_file.read(reinterpret_cast<char *>(&value), sizeof(float));
          facet.at(j)(k) = value;
        }

        add(facet.at(j));
      }

      input_file.ignore(2);
    }
  }

  octree = Octree(*this, 128, 8);

  for (const Facet &facet: facets) {
    octree.add(facet);
  }

  if (update->method_type.compare("tlmpm") == 0) {
    if (domain->boxlo[0] > interval_x.x0)
      domain->boxlo[0] = interval_x.x0;
    if (domain->boxlo[1] > interval_y.x0)
      domain->boxlo[1] = interval_y.x0;
    if (domain->dimension == 3)
      if (domain->boxlo[2] > interval_z.x0)
        domain->boxlo[2] = interval_z.x0;

    if (domain->boxhi[0] < interval_x.x1)
      domain->boxhi[0] = interval_x.x1;
    if (domain->boxhi[1] < interval_y.x1)
      domain->boxhi[1] = interval_y.x1;
    if (domain->dimension == 3)
      if (domain->boxhi[2] < interval_z.x1)
        domain->boxhi[2] = interval_z.x1;
  } else {
    if (domain->boxlo[0] > interval_x.x0)
      interval_x.x0 = domain->boxlo[0];
    if (domain->boxlo[1] > interval_y.x0)
      interval_y.x0 = domain->boxlo[1];
    if (domain->dimension == 3)
      if (domain->boxlo[2] > interval_z.x0)
        interval_z.x0 = domain->boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Stl::inside(double x, double y, double z)
{
  Vector3d p(x, y, z);
  
  double theta = M_PI*rand()/RAND_MAX;
  double phi = 2*M_PI*rand()/RAND_MAX;
  Vector3d direction(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

  return octree.intersections(p, direction).size()%2;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> Stl::limits(){
  vector<double> lim;
  lim.push_back(interval_x.x0);
  lim.push_back(interval_x.x1);
  lim.push_back(interval_y.x0);
  lim.push_back(interval_y.x1);
  lim.push_back(interval_z.x0);
  lim.push_back(interval_z.x1);
  return lim;
}

void Stl::write_restart(ofstream *of) {
  int input_file_name_length = input_file_name.length();
  of->write(reinterpret_cast<const char *>(&input_file_name_length), sizeof(int));
  of->write(input_file_name.c_str(), input_file_name_length);

  of->write(reinterpret_cast<const char *>(&interval_x.x0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&interval_x.x1), sizeof(double));
  of->write(reinterpret_cast<const char *>(&interval_y.x0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&interval_y.x1), sizeof(double));
  of->write(reinterpret_cast<const char *>(&interval_z.x0), sizeof(double));
  of->write(reinterpret_cast<const char *>(&interval_z.x1), sizeof(double));
}

void Stl::read_restart(ifstream *ifr) {
  // cout << "Restart Stl" << endl;
  int input_file_name_length;
  ifr->read(reinterpret_cast<char *>(&input_file_name_length), sizeof(int));
  vector<char> buffer(input_file_name_length + 1);
  ifr->read(buffer.data(), input_file_name_length);
  buffer.back() = '\0';
  input_file_name = string(buffer.data());

  ifr->read(reinterpret_cast<char *>(&interval_x.x0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&interval_x.x1), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&interval_y.x0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&interval_y.x1), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&interval_z.x0), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&interval_z.x1), sizeof(double));
}
