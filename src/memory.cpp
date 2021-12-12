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
#include <string>
#include <string.h>
#include <memory.h>
#include <error.h>

using namespace std;

/* ---------------------------------------------------------------------- */

Memory::Memory(MPM *mpm) : Pointers(mpm) {}

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, string name)
{
  if (nbytes == 0) return nullptr;

#if defined(MPM_MEMALIGN)
  void *ptr;
  int retval = posix_memalign(&ptr, MPM_MEMALIGN, nbytes);
  if (retval) ptr = nullptr;
#else
  void *ptr = malloc(nbytes);
#endif
  if (ptr == nullptr) {
    error->all(FLERR, "Failed to allocate " + to_string(nbytes) + " bytes for array " + name + ".\n");
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, string name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return nullptr;
  }

  ptr = realloc(ptr,nbytes);
  if (ptr == nullptr) {
    error->all(FLERR, "Failed to allocate " + to_string(nbytes) + " bytes for array " + name + ".\n");
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == nullptr) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(string name)
{
  error->all(FLERR, "Cannot create/grow a vector/array of pointers for " + name + ".\n");
}
