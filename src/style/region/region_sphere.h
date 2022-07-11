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

RegionStyle(sphere,Sphere)

#else

#ifndef MPM_REGION_SPHERE_H
#define MPM_REGION_SPHERE_H

#include <region.h>

/*! \ingroup region regionblock region_block

\section Syntax Syntax
\code
region(region-ID, sphere, x, y, z, R)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>x, y, z: coordinates of the sphere center. Only x is required in 1D, and x and y in 2D.</li>
<li>R: sphere radius</li>
</ul>

\section Examples Examples
\code
R = 0.2
region(rSph, sphere, 0, 0, 0, R)
\endcode
Defines a region called 'rSph' delimited by a sphere of radius \f$R=0.2\f$ centered around the origin.

\section Description Description

This command defines a spherical region of space.It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/

class Sphere : public Region {

 public:
  Sphere(class MPM *, vector<string>);
  ~Sphere();
  int inside(float, float, float);
  vector<float> limits();
  void write_restart(ofstream *);
  void read_restart(ifstream *);

 protected:
  float c1, c2, c3, R, RSq;
  float xlo, xhi, ylo, yhi, zlo, zhi;
  string usage[3] = {"Usage in 1D: region(region-ID, sphere, center_x, R)\n",
		     "Usage in 3D: region(region-ID, sphere, center_x, center_y, R)\n",
		     "Usage in 3D: region(region-ID, sphere, center_x, center_y, center_z, R)\n"};
  int Nargs[3] = {4, 5, 6};
};

#endif
#endif
