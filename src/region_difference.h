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

RegionStyle(difference,Difference)

#else

#ifndef MPM_REGION_DIFFERENCE_H
#define MPM_REGION_DIFFERENCE_H

#include "region.h"

/*! \ingroup region regiondifference region_difference

\section Syntax Syntax
\code
region(region-ID, difference, region-1, region-2)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>region-1: name of the region to which region-2 will be subtracted.</li>
<li>region-2: name of the region to be subtracted from region-1.</li>
</ul>

\section Examples Examples
\code
L = 0.5
region(A, block, -L, L, -L, L, -L, L)
region(B, cylinder, x, 0, 0, L, -1, 1)
region(C, difference, A, B)
\endcode
Defines the difference of a rectangle (Block_) called A with a cylinder (Cylinder) called B. This will create a hole in the rectangle A.

\section Description Description

This command defines a difference between two regions of space. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/

class Difference : public Region {

public:
  Difference(class MPM *, vector<string>);
  ~Difference();
  int inside(double, double, double);
  vector<double> limits();
  void write_restart(ofstream *);
  void read_restart(ifstream *);

protected:
  // vector<string> regions;
  vector<int> iregions;
  double xlo, xhi, ylo, yhi, zlo, zhi;

  const string usage =
      "Usage: region(region-ID, difference, region-1, region-2)\n";
  const int Nargs = 4;
};

#endif
#endif
