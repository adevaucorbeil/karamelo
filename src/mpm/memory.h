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
 * -----------------------------------------------------------------------  */

#ifndef MPM_MEMORY_H
#define MPM_MEMORY_H

#include <mpmtype.h>
#include <pointers.h>

/*! This class handles the allocation of all large vectors and arrays.
 * This class is directly herited (aka copied) from LAMMPS.\n\n
 * Create/grow/destroy vecs and multidim arrays with contiguous memory blocks
 * only use with primitive data types, e.g. 1d vec of ints, 2d array of doubles
 * fail() prevents use with pointers, 
 * e.g. 1d vec of int*, due to mismatched destroy
 * avoid use with non-primitive data types to avoid code bloat
 * for these other cases, use smalloc/srealloc/sfree directly
 */
class Memory : protected Pointers {
 public:
  Memory(class MPM *);

  void *smalloc(bigint n, string);           ///< Safe malloc
  void *srealloc(void *, bigint n, string);  ///< Safe reallocation
  void sfree(void *);                        ///< Safe free
  void fail(string);                         ///< Erroneous usage of templated create/grow functions

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE *create(TYPE *&array, int n, string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n;
    array = (TYPE *) smalloc(nbytes,name);
    return array;
  }

  template <typename TYPE>
  TYPE **create(TYPE **&array, int n, string name) 
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   grow or shrink 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE *grow(TYPE *&array, int n, string name)
  {
    if (array == nullptr) return create(array,n,name);

    bigint nbytes = ((bigint) sizeof(TYPE)) * n;
    array = (TYPE *) srealloc(array,nbytes,name);
    return array;
  }

  template <typename TYPE>
  TYPE **grow(TYPE **&array, int n, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   destroy a 1d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy(TYPE *array)
  {sfree(array);}

/* ----------------------------------------------------------------------
   create a 1d array with index from nlo to nhi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE *create1d_offset(TYPE *&array, int nlo, int nhi, string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * (nhi-nlo+1);
    array = (TYPE *) smalloc(nbytes,name);
    array -= nlo;
    return array;
  }

  template <typename TYPE>
  TYPE **create1d_offset(TYPE **&array, int nlo, int nhi, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   destroy a 1d array with index offset
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy1d_offset(TYPE *array, int offset)
  {
    if (array) sfree(&array[offset]);
  }

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE **create(TYPE **&array, int n1, int n2, string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes,name);
    
    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  template <typename TYPE>
  TYPE ***create(TYPE ***&array, int n1, int n2, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE **grow(TYPE **&array, int n1, int n2, string name)
  {
    if (array == nullptr) return create(array,n1,n2,name);
    
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
    TYPE *data = (TYPE *) srealloc(array[0],nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) srealloc(array,nbytes,name);
    
    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }
  
  template <typename TYPE>
  TYPE ***grow(TYPE ***&array, int n1, int n2, string name)
  {fail(name); return nullptr;}
  
/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy(TYPE **array)
  {
    if (array == nullptr) return;
    sfree(array[0]);
    sfree(array);
  }

/* ----------------------------------------------------------------------
   create a 2d array with a ragged 2nd dimension
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE **create_ragged(TYPE **&array, int n1, int *n2, string name)
  {
    bigint n2sum = 0;
    for (int i = 0; i < n1; i++) n2sum += n2[i];
    
    bigint nbytes = ((bigint) sizeof(TYPE)) * n2sum;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes,name);
    
    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2[i];
    }
    return array;
  }
  
  template <typename TYPE>
  TYPE ***create_ragged(TYPE ***&array, int n1, int *n2, string name)
  {fail(name); return nullptr;}
  
/* ----------------------------------------------------------------------
   create a 2d array with 2nd index from n2lo to n2hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE **create2d_offset(TYPE **&array, int n1, int n2lo, int n2hi,
                         string name)
  {
    int n2 = n2hi - n2lo + 1;
    create(array,n1,n2,name);
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    return array;
  }
  
  template <typename TYPE>
  TYPE ***create2d_offset(TYPE ***&array, int n1, int n2lo, int n2hi,
                          string name) {fail(name); return nullptr;}
  
/* ----------------------------------------------------------------------
   destroy a 2d array with 2nd index offset
------------------------------------------------------------------------- */
  
  template <typename TYPE>
  void destroy2d_offset(TYPE **array, int offset)
  {
    if (array == nullptr) return;
    sfree(&array[0][offset]);
    sfree(array);
  }
  
/* ----------------------------------------------------------------------
   create a 3d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ***create(TYPE ***&array, int n1, int n2, int n3, string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
    TYPE **plane = (TYPE **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1;
    array = (TYPE ***) smalloc(nbytes,name);
    
    int i,j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }

  template <typename TYPE>
  TYPE ****create(TYPE ****&array, int n1, int n2, int n3, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, string name)
  {
    if (array == nullptr) return create(array,n1,n2,n3,name);
    
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
    TYPE *data = (TYPE *) srealloc(array[0][0],nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
    TYPE **plane = (TYPE **) srealloc(array[0],nbytes,name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1;
    array = (TYPE ***) srealloc(array,nbytes,name);
    
    int i,j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint) i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m+j] = &data[n];
        n += n3;
      }
    }
    return array;
  }
  
  template <typename TYPE>
  TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   destroy a 3d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy(TYPE ***array)
  {
    if (array == nullptr) return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
  }

/* ----------------------------------------------------------------------
   create a 3d array with 1st index from n1lo to n1hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi,
                          int n2, int n3, string name)
  {
    int n1 = n1hi - n1lo + 1;
    create(array,n1,n2,n3,name);
    array -= n1lo;
    return array;
  }
  
  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi,
                           int n2, int n3, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   free a 3d array with 1st index offset
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy3d_offset(TYPE ***array, int offset)
  {
    if (array) destroy(&array[offset]);
  }

/* ----------------------------------------------------------------------
   create a 3d array with
   1st index from n1lo to n1hi inclusive,
   2nd index from n2lo to n2hi inclusive,
   3rd index from n3lo to n3hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi,
                          int n2lo, int n2hi, int n3lo, int n3hi,
                          string name)
  {
    int n1 = n1hi - n1lo + 1;
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    create(array,n1,n2,n3,name);
    
    bigint m = ((bigint) n1) * n2;
    for (bigint i = 0; i < m; i++) array[0][i] -= n3lo;
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    array -= n1lo;
    return array;
  }
  
  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****&array, int n1lo, int n1hi,
                           int n2lo, int n2hi, int n3lo, int n3hi,
                           string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   free a 3d array with all 3 indices offset
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy3d_offset(TYPE ***array,
                        int n1_offset, int n2_offset, int n3_offset)
  {
    if (array == nullptr) return;
    sfree(&array[n1_offset][n2_offset][n3_offset]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(&array[n1_offset]);
  }
  
/* ----------------------------------------------------------------------
   create a 4d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4,
                  string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1*n2*n3;
    TYPE **cube = (TYPE **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1*n2;
    TYPE ***plane = (TYPE ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE ***)) * n1;
    array = (TYPE ****) smalloc(nbytes,name);

    int i,j,k;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 + j*n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 + j*n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }
  
  template <typename TYPE>
  TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4,
                   string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   destroy a 4d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy(TYPE ****array)
  {
    if (array == nullptr) return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
  }

/* ----------------------------------------------------------------------
   create a 4d array with indices
   2nd index from n2lo to n2hi inclusive
   3rd index from n3lo to n3hi inclusive
   4th index from n4lo to n4hi inclusive
   cannot grow it
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE ****create4d_offset(TYPE ****&array, int n1, int n2lo, int n2hi,
                           int n3lo, int n3hi, int n4lo, int n4hi,
                           string name)
  {
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    int n4 = n4hi - n4lo + 1;
    create(array,n1,n2,n3,n4,name);
    
    bigint m = ((bigint) n1) * n2 * n3;
    for (bigint i = 0; i < m; i++) array[0][0][i] -= n4lo;
    m = ((bigint) n1) * n2;
    for (bigint i = 0; i < m; i++) array [0][i] -= n3lo;
    for (int i = 0; i < n1; i++) array[i] -= n2lo;
    return array;
  }
  
  template <typename TYPE>
  TYPE ****create4d_offset(TYPE *****&array, int n1, int n2lo, int n2hi,
                           int n3lo, int n3hi, int n4lo, int n4hi,
                           string name)
  {fail(name); return nullptr;}
  
/* ----------------------------------------------------------------------
   free a 4d array with indices 2,3, and 4 offset
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy4d_offset(TYPE ****array, 
                        int n2_offset, int n3_offset, int n4_offset)
  {
    if (array == nullptr) return;
    sfree(&array[0][n2_offset][n3_offset][n4_offset]);
    sfree(&array[0][n2_offset][n3_offset]);
    sfree(&array[0][n2_offset]);
    sfree(array);
  }

/* ----------------------------------------------------------------------
   create a 5d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4,
                   int n5, string name)
  {
    bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4*n5;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE *)) * n1*n2*n3*n4;
    TYPE **level4 = (TYPE **) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE **)) * n1*n2*n3;
    TYPE ***level3 = (TYPE ***) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE ***)) * n1*n2;
    TYPE ****level2 = (TYPE ****) smalloc(nbytes,name);
    nbytes = ((bigint) sizeof(TYPE ****)) * n1;
    array = (TYPE *****) smalloc(nbytes,name);
    
    int i,j,k,l;
    bigint m1,m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint) i) * n2;
      array[i] = &level2[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint) i) * n2 + j;
        m2 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3;
        level2[m1] = &level3[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint) i) * n2*n3 +  ((bigint) j) * n3 + k;
          m2 = ((bigint) i) * n2*n3*n4 +
            ((bigint) j) * n3*n4 + ((bigint) k) * n4;
          level3[m1] = &level4[m2];
          for (l = 0; l < n4; l++) {
            m1 = ((bigint) i) * n2*n3*n4 +
              ((bigint) j) * n3*n4 + ((bigint) k) * n4 + l;
            level4[m1] = &data[n];
            n += n5;
          }
        }
      }
    }
    return array;
  }
  
  template <typename TYPE>
  TYPE ******create(TYPE ******&array, int n1, int n2, int n3, int n4,
                    int n5, string name)
  {fail(name); return nullptr;}

/* ----------------------------------------------------------------------
   destroy a 5d array
------------------------------------------------------------------------- */

  template <typename TYPE>
  void destroy(TYPE *****array)
  {
    if (array == nullptr) return;
    sfree(array[0][0][0][0]);
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
  }
  
/* ----------------------------------------------------------------------
   memory usage of arrays, including pointers
------------------------------------------------------------------------- */

  template <typename TYPE>
  bigint usage(TYPE *array, int n)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n;
    return bytes;
  }
  
  template <typename TYPE>
  bigint usage(TYPE **array, int n1, int n2)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2;
    bytes += ((bigint) sizeof(TYPE *)) * n1;
    return bytes;
  }
  
  template <typename TYPE>
  bigint usage(TYPE ***array, int n1, int n2, int n3)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
    bytes += ((bigint) sizeof(TYPE *)) * n1*n2;
    bytes += ((bigint) sizeof(TYPE **)) * n1;
    return bytes;
  }
  
  template <typename TYPE>
  bigint usage(TYPE ****array, int n1, int n2, int n3, int n4)
  {
    (void) array;
    bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2*n3*n4;
    bytes += ((bigint) sizeof(TYPE *)) * n1*n2*n3;
    bytes += ((bigint) sizeof(TYPE **)) * n1*n2;
    bytes += ((bigint) sizeof(TYPE ***)) * n1;
    return bytes;
  }
};


#endif
