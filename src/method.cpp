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

#include "method.h"

using namespace std;


Method::Method(MPM *mpm) : Pointers(mpm)
{
  is_TL = false;
  is_CPDI = false;
  ge = false;
}

Method::~Method()
{
}
