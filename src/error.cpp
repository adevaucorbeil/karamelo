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
#include "string.h"
#include "universe.h"
#include "input.h"
#include "error.h"

using namespace std;


Error::Error(MPM *mpm) : Pointers(mpm) {
}


/*! Called by all procs in universe
 * close all output, screen, and log files in universe
 * no abort, so insure all procs in universe call, else will hang
 */


void Error::all(const char *file, int line, const string str)
{
  MPI_Barrier(universe->uworld);

  int me;
  const char *lastcmd = (const char*)"(unknown)";

  MPI_Comm_rank(universe->uworld, &me);

  if (universe->me == 0) {
    cout << "Error at line " << input->line_number << ": " << str
         << " raised at (" << file << "," << line << ")\n";
    cout << "Last command: " << input->line << endl;
  }

  MPI_Finalize();
  exit(1);
}

/*! Called by one proc in universe
 *  Forces abort of entire universe if any proc in universe calls
 */

void Error::one(const char *file, int line, const string str)
{
  int me;
  const char *lastcmd = (const char*)"(unknown)";

  MPI_Comm_rank(universe->uworld, &me);

  cout << "Error at line " << input->line_number << ": " << str
       << " raised at (" << file << "," << line << ")\n";
  cout << "Last command: " << input->line << endl;

  MPI_Abort(universe->uworld, 1);
}

void Error::done(int status)
{
  MPI_Barrier(universe->uworld);

  MPI_Finalize();
  exit(status);
}
