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

#ifndef LMP_GROUP_H
#define LMP_GROUP_H

#include "pointers.h"
#include <string>
#include <vector>

using namespace std;

class Group : protected Pointers
{
public:
  int ngroup;       ///< Number of defined groups
  string *names;    ///< Name of each group
  int *bitmask;     ///< One-bit mask for each group
  //int *inversemask; ///< Inverse mask for each group
  string *pon;      ///< Group of particles if pon == "particles", or nodes if pon = "nodes"
  int *solid;       ///< Solids corresponding to each group
  int *region;      ///< Index of the region corresponding to each group

  Group(class MPM *);
  virtual ~Group();

  void assign(vector<string>); ///< Assign atoms to a new or existing group
  int find(string);            ///< Return group index
  int find_unused();           ///< Return index of first available group

  double xcm(int, int);        ///< Determine the centre of mass of a group
  double internal_force(
      int,
      int);                    ///< Determine the resulting internal force applied onto the group
  double external_force(
      int,
      int);                    ///< Determine the resulting internal force applied onto the group

  void write_restart(ofstream *);
  void read_restart(ifstream *);  
};

#endif

/*! \defgroup group group

\section Syntax Syntax
\code
group(group-ID, nodes, region, region-ID, solid, solid-ID)
group(group-ID, particles, region, region-ID, solid, solid-ID)
\endcode

Select a group of nodes or particles from a given solid (solid-ID) that lie within a given region (region-ID) at the start of the simulation.


\section Examples Examples
\code
hLx = 1
solid(solid1, region, box, ppc, mat1, cellsize, Tr)
region(sidexmax, block, hLx - 0.6 * cellsize, INF, INF, INF)
group(sidexmax_n, nodes, region, sidexmax, solid, solid1)
\endcode
Defines a group of nodes named 'sidexmax_n' from the background grid of solid 'solid1' which lie in the region 'sidexmax' as defined just above. If the method is TLMPM, TLCPDI or TLCPDI2, each solid has its own background grid, this is why it is important to specify the solid. If using updated Lagrangian, i.e. ULMPM, ULCPDI, and ULCPDI2, there is only one background grid, so any solid specify would give the same result. However, in these cases, the last two parameters are not optional.

\section Description Description

This command creates a group of nodes or particles lying within a region at the beginning of the simulation. If these particles move out of the region during the simulation, they stay in the group. If new particles move into the region during the simulation, they will not be part of the group.

*/
