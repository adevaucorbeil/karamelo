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

#include <iostream>
#include "damage.h"
#include "error.h"

using namespace std;


Damage::Damage(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new Damage with ID: " << args[0] << endl;
  id = args[0];
}

Damage::~Damage()
{
}


void Damage::init()
{
}

void Damage::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In Damage::options()" << endl;

  if (args->end() < it) error->all(FLERR, "Error: not enough arguments\n");

  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

