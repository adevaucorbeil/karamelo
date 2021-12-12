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

#include <eos.h>
#include <error.h>
#include <universe.h>
#include <iostream>

using namespace std;


EOS::EOS(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  if (universe->me == 0) {
    cout << "Creating new EOS with ID: " << args[0] << endl;
  }
  id = args[0];
  style = args[1];
}

EOS::~EOS()
{
}


void EOS::init()
{
}

void EOS::options(vector<string> *args, vector<string>::iterator it)
{
  if (universe->me == 0) {
    cout << "In EOS::options()" << endl;
  }
  if (args->end() < it) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }
  if (args->end() > it) {
    if (universe->me == 0) {
      cout << "Ignoring optional arguments: ";
      for (it; it != args->end(); ++it){
	cout << *it << "\t";
      }
      cout << endl;
    }
  }
}

