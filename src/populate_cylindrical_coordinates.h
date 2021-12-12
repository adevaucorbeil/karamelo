/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(populate_cylindrical_coordinates,PopulateCylindricalCoordinates)

#else

#ifndef MPM_POPULATE_CYLINDRICAL_COORDINATE_H
#define MPM_POPULATE_CYLINDRICAL_COORDINATE_H

#include <pointers.h>
#include <var.h>

class PopulateCylindricalCoordinates : protected Pointers {
 public:
  PopulateCylindricalCoordinates(class MPM *);
  class Var command(vector<string>);

private:
  const map<int, string> usage = {
      {
          2,
          "Usage: populate_cylindrical_coordinates(solid-ID, region-ID, c1, c2, R, Ndr, Ndtheta, T0)\n",
      },
      {
          3,
          "Usage: populate_cylindrical_coordinates(solid-ID, region-ID, "
          "axis, c1, c2, R, lo, hi, Ndx, Ndr, Ndtheta, T0)\n",
      }};
  const map<int, int> Nargs = {{2, 8}, {3, 12}};
};

#endif
#endif
