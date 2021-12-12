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

#ifndef MPM_ERROR_H
#define MPM_ERROR_H

#include <pointers.h>

/*! This class handles all the errors that can occur during the execution of the code.
 *  Terminates the process and prints all the corresponding error messages.
 */
class Error : protected Pointers {
 public:
  Error(class MPM *);

  void all(const char *, int, const string);      ///< Is used when all CPUs exit with an error.
  void one(const char *, int, const string);      ///< Is used when one CPU exits with an error.
  void done(int = 0);                             ///< Exit with status
};

#endif
