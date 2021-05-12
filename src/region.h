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

/* \page region
 */

#ifndef MPM_REGION_H
#define MPM_REGION_H

#include "pointers.h"
#include <vector>

/*! Parent class of all the different kinds of regions that can be used.
 *
 * Stores the region id as well as functions that returns the regions' 
 * limits and on that checks if a point lies within the region's limits.
 */
class Region : protected Pointers {
 public:
  string id;                        ///< Region identification string
  string style;                     ///< Region style
  int interior;                     ///< 1 for interior, 0 for exterior

  Region(class MPM *, vector<string>);
  virtual ~Region();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // called by other classes to check point versus region

  int match(double, double, double);

  // implemented by each region

  virtual vector<double> limits() {return vector<double>();};
  virtual int inside(double, double, double) = 0;
  virtual void write_restart(ofstream*) = 0;
  virtual void read_restart(ifstream*) = 0;
  //protected:
};

#endif

/*! \defgroup region region

\section Syntax Syntax
\code
region(region-ID, region_type, region_specific_arguments)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>region-type: block, cylinder, sphere, ...</li>
<li>region_specific_arguments: list of arguments specific to the region type used. </li>
</ul>

\section Examples Examples
\code
R = 0.2
region(rSph, sphere, 0, 0, 0, R)
\endcode
Defines a region called 'rSph' delimited by a sphere of radius \f$R=0.2\f$ centered around the origin.

\code
region(A, block, 0, L, 0, 2*L)
\endcode
Defines a 2D rectangular region called 'A'. A point of coordinates (x, y, z) lies within the region if \f$0 \leq x \leq L\f$ and \f$0 \leq y \leq 2L\f$.

\section Description Description

This command defines a geometric region of space. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Region_Subclasses Region types

To access the different region type supported, and how to use the region() command with them, refer to the corresponding class.
*/
