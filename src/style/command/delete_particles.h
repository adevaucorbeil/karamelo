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

CommandStyle(delete_particles,DeleteParticles)

#else

#ifndef MPM_DELETE_PARTICLES_H
#define MPM_DELETE_PARTICLES_H

#include <vector>
#include <pointers.h>
#include <deque>

class DeleteParticles : protected Pointers {
 public:
  DeleteParticles(class MPM *);
  class Var command(vector<string>);


  void delete_region(vector<string>, int);
};

#endif
#endif
