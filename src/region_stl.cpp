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

#define BIG 1.0e20

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

  fstream input_file(input_file_name.c_str(), fstream::in);

  if (!input_file) {
    error->all(FLERR, "Error: input file could not be opened.\n");
  }

  xlo =  BIG;
  xhi = -BIG;
  ylo =  BIG;
  yhi = -BIG;
  zlo =  BIG;
  zhi = -BIG;

  string token;
  
  cin >> token;

  if (token != "solid") {
    error->all(FLERR, "Error: stl file not valid.\n");
  }

  // note: name might be optional?
  cin >> name;

  while (true) {
    cin >> token;

    if (token == "endsolid") {
      break;
    }

    if (token != "facet") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }

    cin >> token;

    if (token != "normal") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }

    double normal[3];

    for (int i = 0; i < 3; i++) {
      cin >> normal[i];
    }

    cin >> token;

    if (token != "outer") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }

    cin >> token;

    if (token != "loop") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }

    for (int i = 0; i < 3; i++) {
      cin >> token;

      if (token != "vertex") {
        error->all(FLERR, "Error: stl file not valid.\n");
      }

      double vertex[3];

      for (int j = 0; j < 3; j++) {
        cin >> vertex[j];
      }

      xlo = min(xlo, vertex[0]);
      xhi = max(xhi, vertex[0]);
      ylo = min(xlo, vertex[1]);
      yhi = max(xhi, vertex[1]);
      zlo = min(xlo, vertex[2]);
      zhi = max(xhi, vertex[2]);
    }

    cin >> token;

    if (token != "endloop") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }

    cin >> token;

    if (token != "endfacet") {
      error->all(FLERR, "Error: stl file not valid.\n");
    }
  }


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
  return 0;
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
