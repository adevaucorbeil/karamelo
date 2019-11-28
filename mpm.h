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

#ifndef MPM_MPM_H
#define MPM_MPM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <mpi.h>

using namespace std;

class MPM {

 public:

  class Memory *memory;          // memory allocation functions
  class Error *error;            // error handling
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
  class Output *output;          // thermo/dump/restart

  class Domain *domain;          // simulation box
  class Material *material;      // material
  class Update *update;          //
  class Modify *modify;          // fixes and computes
  class Group *group;            // groups of particles

  MPI_Comm world;                // MPI communicator
  double initclock;              // wall clock at instantiation

  filebuf infile;                // infile
  filebuf logfile;               // logfile
  ofstream *wlogfile;

  MPM(int, char **, MPI_Comm);
  ~MPM();
  void init();

private:
  void help();
  MPM() {};                   // prohibit using the default constructor
  MPM(const MPM &) {};        // prohibit using the copy constructor
};

#endif
