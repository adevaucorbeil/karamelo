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

RegionStyle(stl,Stl)

#else

#ifndef MPM_REGION_STL_H
#define MPM_REGION_STL_H

#include "region.h"

/*! \ingroup region regionblock region_block

\section Syntax Syntax
\code
region(region-ID, stl, input_file.stl)
\endcode

<ul>
<li>region-ID: name of the region to be created.</li>
<li>input_file.stl: path to stl file.</li>
</ul>

\section Examples Examples

\section Description Description

This command defines a region of space described by an stl file. It is usually used by groups to find what nodes and/or particles lie within a specific region of space. A region does not select nodes or particles.

\section Class Class description
*/

class Stl : public Region {
 public:
  Stl(class MPM *, vector<string>);
  int inside(double, double, double);
  vector<double> limits();
  void write_restart(ofstream *);
  void read_restart(ifstream *);

 protected:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  string input_file_name;
  string name;
  string usage[3] = {"Usage: region(region-ID, stl, input_file.stl)\n"};
};

#endif
#endif
