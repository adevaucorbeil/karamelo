#include <iostream>
#include "string.h"
#include "universe.h"
#include "input.h"
#include "error.h"

using namespace std;


Error::Error(MPM *mpm) : Pointers(mpm) {
}



void Error::all(const char *file, int line, const string str)
{
  MPI_Barrier(universe->uworld);

  int me;
  const char *lastcmd = (const char*)"(unknown)";

  MPI_Comm_rank(universe->uworld,&me);

  // if (me == 0) {
  //   if (input && input->line) lastcmd = input->line;
  //   if (screen) fprintf(screen,"ERROR: %s (%s:%d)\n"
  //                       "Last command: %s\n",
  //                       str,file,line,lastcmd);
  //   if (logfile) fprintf(logfile,"ERROR: %s (%s:%d)\n"
  //                        "Last command: %s\n",
  //                        str,file,line,lastcmd);
  // }
  cout << "Error: " << str << "(" << file << "," << line << ")\n";
  cout << "Last command: " << input->line << endl;
  MPI_Finalize();
  exit(1);
}
