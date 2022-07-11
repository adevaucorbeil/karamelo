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

#include <region_cylinder.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <math_special.h>
#include <method.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <iostream>

using namespace std;
using namespace MathSpecial;

#define BIG 1.0e20

Cylinder::Cylinder(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    c1 = c2 = R = RSq = lo = hi = 0;
    axis = 'x';
    xlo = xhi = ylo = yhi = zlo = zhi = 0;
    return;
  }

  if (universe->me == 0)
    cout << "Initiate Cylinder" << endl;

  if (domain->dimension == 3) {
    if (args.size()<8) {
      error->all(FLERR, "Error: not enough arguments.\n");
    }

    options(&args, args.begin()+8);
  } else {
    if (args.size()<5) {
      error->all(FLERR, "Error: not enough arguments.\n");
    }    

    options(&args, args.begin()+5);
  }

  if (domain->dimension == 3) {
    if (args[2] == "x") {
      axis = 'x';
    } else if (args[2] == "y") {
      axis = 'y';
    } else if (args[2] == "z") {
      axis = 'z';
    } else {
      error->all(FLERR, "Error: region cylinder axis not understood, expect x, y, or z, received " + args[2]+".\n");
    }

    c1 = input->parsev(args[3]);
    c2 = input->parsev(args[4]);
    R = input->parsev(args[5]);

    if (args[6] == "INF" || args[6] == "EDGE") {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      lo = -BIG;
    } else lo = input->parsev(args[6]);

    if (args[7] == "INF" || args[7] == "EDGE") {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      hi = BIG;
    } else hi = input->parsev(args[7]);

    if (universe->me == 0)
      cout << "lo hi = " << lo << "\t" << hi << endl;
  } else {
    axis = 'z';

    c1 = input->parsev(args[2]);
    c2 = input->parsev(args[3]);
    R = input->parsev(args[4]);
    lo = hi = 0;
  }

  RSq = R*R;
  if (universe->me == 0)
    cout << "axis, c1, c2, R = " << axis << "\t" << c1 << "\t" << c2 << "\t" << R << endl;

  // error check

  if (lo > hi) {
    error->all(FLERR, "Illegal region cylinder command: low is higher than high.\n");
  }

  if (axis=='x') {
    xlo = lo;
    xhi = hi;
    ylo = c1 - R;
    yhi = c1 + R;
    zlo = c2 - R;
    zhi = c2 + R;
  } else if (axis=='y') {
    xlo = c1 - R;
    xhi = c1 + R;
    ylo = lo;
    yhi = hi;
    zlo = c2 - R;
    zhi = c2 + R;
  } else if (axis=='z') {
    xlo = c1 - R;
    xhi = c1 + R;
    ylo = c2 - R;
    yhi = c2 + R;
    zlo = lo;
    zhi = hi;
  }
}


Cylinder::~Cylinder()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Cylinder::inside(float x, float y, float z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  float dSq;
  if (axis=='x') {
    dSq = square(y - c1) + square(z - c2);
    if (x >= lo && x <= hi && dSq <= RSq) {
      return 1;
    }
  } else if (axis=='y') {
    dSq = square(x - c1) + square(z - c2);
    if (y >= lo && y <= hi && dSq <= RSq) {
      return 1;
    }
  } else if (axis=='z') {
    dSq = square(x - c1) + square(y - c2);
    if (z >= lo && z <= hi && dSq <= RSq) {
      return 1;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<float> Cylinder::limits(){
  vector<float> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}

void Cylinder::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&c1), sizeof(float));
  of->write(reinterpret_cast<const char *>(&c2), sizeof(float));
  of->write(reinterpret_cast<const char *>(&R), sizeof(float));
  of->write(reinterpret_cast<const char *>(&lo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&hi), sizeof(float));
  of->write(&axis, sizeof(char));
  of->write(reinterpret_cast<const char *>(&xlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&xhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&ylo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&yhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zhi), sizeof(float));
}

void Cylinder::read_restart(ifstream *ifr)
{
  // cout << "Restart Cylinder" << endl;
  ifr->read(reinterpret_cast<char *>(&c1), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&c2), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&R), sizeof(float));
  RSq = R * R;
  ifr->read(reinterpret_cast<char *>(&lo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&hi), sizeof(float));
  ifr->read(&axis, sizeof(char));
  ifr->read(reinterpret_cast<char *>(&xlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&xhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&ylo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&yhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zhi), sizeof(float));
}
