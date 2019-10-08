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
#include "dump.h"
#include "error.h"

using namespace std;


Dump::Dump(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new dump with ID: " << args[0] << endl;
  id = args[0];

  if (args[1].compare("all")!=0) {
    error->all(FLERR, "Error: groups are not yet supported for dumps. Should use all.\n");
  }

  style = args[2];
  
  filename = args[4];
}

Dump::~Dump()
{
}


void Dump::options(vector<string> *args, vector<string>::iterator it)
{
  cout << "In dump::options()" << endl;
  if (args->end() < it) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it) {
    cout << "Ignoring optional arguments: ";
    for (it; it != args->end(); ++it){
      cout << *it << "\t";
    }
    cout << endl;
  }
}

