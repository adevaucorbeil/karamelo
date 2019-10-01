/* -*- c++ -*- ----------------------------------------------------------*/

#include <mpi.h>
#include <string>

using namespace std;

inline void MPI_string_bcast(string &data, MPI_Datatype datatype, int root, MPI_Comm communicator) {
  int n, me;

  MPI_Comm_rank(communicator,&me);

  n = data.size();
  MPI_Bcast(&n,1,MPI_INT,root,communicator);


  char * buf = new char[n + 1];
  copy(data.begin(), data.end(), buf);
  buf[data.size()] = '\0'; // don't forget the terminating 0

  MPI_Bcast(buf, n+1, MPI_CHAR, root, communicator);
  if (me != root) data = string(buf);
  
  delete[] buf;
};
