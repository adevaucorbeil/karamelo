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

#include "damage.h"
#include "error.h"
#include "universe.h"
#include <iostream>

using namespace std;


Damage::Damage(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  if (universe->me == 0)
    cout << "Creating new Damage with ID: " << args[0] << endl;
  id = args[0];
  style = args[1];
}

Damage::~Damage()
{
}


void Damage::init()
{
}

void Damage::options(vector<string> *args, vector<string>::iterator it)
{
  // cout << "In Damage::options()" << endl;

  if (args->end() < it) error->all(FLERR, "Error: not enough arguments\n");

  if (args->end() > it && universe->me == 0) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

