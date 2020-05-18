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

RegionStyle(block,Block_)

#else

#ifndef MPM_REGION_BLOCK_H
#define MPM_REGION_BLOCK_H

#include "region.h"

/*! \ingroup region regionblock region_block

\section Syntax Syntax
\code
region(region-ID, block, xmin, xmax, ymin, ymax, zmin, zmax)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>xmin, xmax: lower and upper boundaries along the x axis.</li>
<li>ymin, ymax: lower and upper boundaries along the y axis. Optional if dimension is 1D.</li>
<li>zmin, zmax: lower and upper boundaries along the z axis. Optional if dimension is 1D or 2D.</li>
</ul>

\section Examples Examples
\code
region(A, block, 0, L, 0, 2*L)
\endcode
Defines a 2D rectangular region called 'A'. A point of coordinates (x, y, z) lies within the region if \f$0 \leq x \leq L\f$ and \f$0 \leq y \leq 2L\f$.

\section Description Description

This command defines a rectangle region of space. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/


class Block_ : public Region {

 public:
  Block_(class MPM *, vector<string>);
  ~Block_();
  int inside(double, double, double);
  vector<double> limits();

 protected:
  double xlo,xhi,ylo,yhi,zlo,zhi;
};

#endif
#endif
