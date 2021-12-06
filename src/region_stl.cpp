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

#include "region_stl.h"
#include "error.h"
#include "domain.h"
#include "universe.h"
#include "update.h"

#include <sstream>

#define BIG 1.0e20
#define EPSILON 0.0000001

bool Facet::intersects(const Vector3d &origin, const Vector3d &direction) const {
  Vector3d edge0 = at(1) - at(0);
  Vector3d edge1 = at(2) - at(0);
  Vector3d h = direction.cross(edge1);
  double a = edge0.dot(h);
  if (abs(a) < EPSILON)
    return false;    // This ray is parallel to this triangle.
  double f = 1/a;
  Vector3d s = origin - at(0);
  double u = f*s.dot(h);
  if (u < 0.0 || u > 1.0)
    return false;
  Vector3d q = s.cross(edge0);
  double v = f*direction.dot(q);
  if (v < 0.0 || u + v > 1.0)
    return false;
  // At this stage we can compute t to find out where the intersection point is on the line.
  double t = f*edge1.dot(q);
  if (t > EPSILON) // ray intersection
  {
    return true;
  }
  // This means that there is a line intersection but not a ray intersection.
  return false;
}

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

  fstream input_file(input_file_name.c_str(), ios::in | ios::binary);

  if (!input_file) {
    error->all(FLERR, "Error: input file could not be opened.\n");
  }

  xlo =  BIG;
  xhi = -BIG;
  ylo =  BIG;
  yhi = -BIG;
  zlo =  BIG;
  zhi = -BIG;

  char buffer[5];
  input_file.read(buffer, 5);

  // whitespace?

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

        xlo = min(xlo, facet.at(i)(0));
        xhi = max(xhi, facet.at(i)(0));
        ylo = min(xlo, facet.at(i)(1));
        yhi = max(xhi, facet.at(i)(1));
        zlo = min(xlo, facet.at(i)(2));
        zhi = max(xhi, facet.at(i)(2));
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

        xlo = min(xlo, facet.at(j)(0));
        xhi = max(xhi, facet.at(j)(0));
        ylo = min(xlo, facet.at(j)(1));
        yhi = max(xhi, facet.at(j)(1));
        zlo = min(xlo, facet.at(j)(2));
        zhi = max(xhi, facet.at(j)(2));
      }

      input_file.ignore(2);
    }
  }

  cout << xlo << ", "
       << xhi << ", "
       << ylo << ", "
       << yhi << ", "
       << zlo << ", "
       << zhi << endl;

  if (update->method_type.compare("tlmpm") == 0) {
    if (domain->boxlo[0] > xlo)
      domain->boxlo[0] = xlo;
    if (domain->boxlo[1] > ylo)
      domain->boxlo[1] = ylo;
    if (domain->dimension == 3)
      if (domain->boxlo[2] > zlo)
        domain->boxlo[2] = zlo;

    if (domain->boxhi[0] < xhi)
      domain->boxhi[0] = xhi;
    if (domain->boxhi[1] < yhi)
      domain->boxhi[1] = yhi;
    if (domain->dimension == 3)
      if (domain->boxhi[2] < zhi)
        domain->boxhi[2] = zhi;
  } else {
    if (domain->boxlo[0] > xlo)
      xlo = domain->boxlo[0];
    if (domain->boxlo[1] > ylo)
      ylo = domain->boxlo[1];
    if (domain->dimension == 3)
      if (domain->boxlo[2] > zlo)
        zlo = domain->boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Stl::inside(double x, double y, double z)
{
  //cout << "testing " << x << ", " << y << ", " << z << endl;
  int intersections = 0;

  Vector3d p(x, y, z);

  for (const Facet &facet: facets) {
    double theta = M_PI*rand()/RAND_MAX;
    double phi = 2*M_PI*rand()/RAND_MAX;
    Vector3d direction(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    if (facet.intersects(p, Vector3d(0, 0, 1))) {
      intersections++;
    }
  }

  return intersections%2;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> Stl::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}

void Stl::write_restart(ofstream *of) {
  int input_file_name_length = input_file_name.length();
  of->write(reinterpret_cast<const char *>(&input_file_name_length), sizeof(int));
  of->write(input_file_name.c_str(), input_file_name_length);

  of->write(reinterpret_cast<const char *>(&xlo), sizeof(double));
  of->write(reinterpret_cast<const char *>(&xhi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&ylo), sizeof(double));
  of->write(reinterpret_cast<const char *>(&yhi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&zlo), sizeof(double));
  of->write(reinterpret_cast<const char *>(&zhi), sizeof(double));
}



void Stl::read_restart(ifstream *ifr) {
  // cout << "Restart Stl" << endl;
  int input_file_name_length;
  ifr->read(reinterpret_cast<char *>(&input_file_name_length), sizeof(int));
  vector<char> buffer(input_file_name_length + 1);
  ifr->read(buffer.data(), input_file_name_length);
  buffer.back() = '\0';
  input_file_name = string(buffer.data());

  ifr->read(reinterpret_cast<char *>(&xlo), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&xhi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&ylo), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&yhi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&zlo), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&zhi), sizeof(double));
}
