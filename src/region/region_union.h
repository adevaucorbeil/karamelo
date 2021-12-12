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


#ifdef REGION_CLASS

RegionStyle(union,Union)

#else

#ifndef MPM_REGION_UNION_H
#define MPM_REGION_UNION_H

#include <region.h>

/*! \ingroup region regionunion region_union

\section Syntax Syntax
\code
region(region-ID, union, region-1, region-2, ..., region-N)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>region-1, region-2, ..., region-N: name of the region to be added.</li>
</ul>

\section Examples Examples
\code
L = 0.5
region(A, block, -L, L, -L, L, -L, L)
region(B, cylinder, x, 0, 0, L, -1, 1)
region(C, union, A, B)
\endcode
Defines the union of a rectangle (Block_) called A and a cylinder (Cylinder) called B.

\section Description Description

This command defines a union between different regions of space. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/


class Union : public Region {

 public:
  Union(class MPM *, vector<string>);
  ~Union();
  int inside(double, double, double);
  vector<double> limits();
  void write_restart(ofstream *);
  void read_restart(ifstream *);

 protected:
  //vector<string> regions;
  vector<int> iregions;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  
  const string usage = "Usage: region(region-ID, union, region-1, region-2, ..., region-N)\n";
  const int Nargs = 4;
};

#endif
#endif
