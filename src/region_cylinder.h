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

RegionStyle(cylinder,Cylinder)

#else

#ifndef MPM_REGION_CYLINDER_H
#define MPM_REGION_CYLINDER_H

#include "region.h"

/*! \ingroup region regioncylinder region_cylinder

\section Syntax Syntax
\code
region(region-ID, cylinder, axis, c1, c2, R, lo, hi)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>axis: direction of the cylinder's axis: 'x', 'y' or 'z'.</li>
<li>c1, c2: coordinates of the cylinder's axis in the other direction.</li>
<li>R: radius of the base of the cylinder.</li>
<li>lo, hi: lower and upper limits of the cylinder in the direction of the axis.</li>
</ul>

\section Examples Examples
\code
region(rCyl, cylinder, x, 0, 0, 10, -1, 1)
\endcode
Defines a cylinder called 'rCyl' of radius 10 with axis x passing through the origin, i.e. (0,0) and spanning between x=-1 to x=+1.

\section Description Description

This command defines a cylindrical region of space. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/
class Cylinder : public Region {

 public:
  Cylinder(class MPM *, vector<string>);
  ~Cylinder();
  int inside(double, double, double);
  vector<double> limits();

 protected:
  double c1, c2, R, RSq, lo, hi;
  char axis;
  double xlo, xhi, ylo, yhi, zlo, zhi;
};

#endif
#endif
