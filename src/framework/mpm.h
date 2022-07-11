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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

#include <mpi.h>

typedef int64_t bigint;
typedef int64_t tagint;

#define MPI_MPM_TAGINT MPI_INT64_T

using namespace std;

/*! Main class creating all the others.*/

class Memory;
class Error;
class Universe;
class Input;
class Output;

class Domain;
class Material;
class Update;
class Modify;
class Group;

class MPM {

 public:

  Memory *memory;          ///< memory allocation functions
  Error *error;            ///< error handling
  Universe *universe;      ///< universe of processors
  Input *input;            ///< input script processing
  Output *output;          ///< logs/dump/restart

  Domain *domain;          ///< simulation box
  Material *material;      ///< material
  Update *update;          ///< pointer to update 
  Modify *modify;          ///< fixes and computes
  Group *group;            ///< groups of particles

  MPI_Comm world;                ///< MPI communicator
  float initclock;              ///< wall clock at instantiation

  filebuf infile;                ///< input file
  //filebuf logfile;               ///< logfile
  ofstream *wlogfile;            ///< log file

  MPM(int, char **, MPI_Comm);   ///< Constructor
  ~MPM();                        ///< Destructor
  void init();

private:
  void help();
  MPM() {};                   ///< prohibit using the default constructor
  MPM(const MPM &) {};        ///< prohibit using the copy constructor
};

#endif
