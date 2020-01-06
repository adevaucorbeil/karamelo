/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMMAND_CLASS

CommandStyle(xcm,CentreOfMass)

#else

#ifndef LMP_CENTRE_OF_MASS_H
#define LMP_CENTRE_OF_MASS_H

#include "pointers.h"

class CentreOfMass : protected Pointers {
 public:
  CentreOfMass(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
