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


#ifndef MPM_POINTERS_H
#define MPM_POINTERS_H

#include "mpm.h"
#include "mpmtype.h"

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

using namespace std;

/*! Father of all classes.*/
class Pointers {
 public:
 Pointers(MPM *ptr) :
   mpm(ptr),
   memory(ptr->memory),
   error(ptr->error),
   universe(ptr->universe),
   input(ptr->input),
   output(ptr->output),
   infile(&ptr->infile),
   //logfile(&ptr->logfile),
   wlogfile(ptr->wlogfile),
   domain(ptr->domain),
   material(ptr->material),
   update(ptr->update),
   modify(ptr->modify),
   group(ptr->group) {}
  virtual ~Pointers() {}
 protected:
  MPM *mpm;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;
  Output *&output;
  
  Domain *&domain;
  Material *&material;
  Update *&update;
  Modify *&modify;
  Group *&group;

  filebuf *infile;
  filebuf *logfile;
  ofstream *&wlogfile;
};

  
#endif
