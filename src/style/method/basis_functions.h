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

#ifndef MPM_BASIS_FUNCTIONS_H
#define MPM_BASIS_FUNCTIONS_H

#include <Kokkos_Core.hpp>

namespace Kokkos
{
  namespace Experimental
  {
    KOKKOS_INLINE_FUNCTION float
    copysign(float mag, float sign)
    {
      return sign == 0? 0:
             sign >  0?  Kokkos::Experimental::abs(mag):
                        -Kokkos::Experimental::abs(mag);
    }
  }
}

namespace BasisFunction
{ 
  KOKKOS_INLINE_FUNCTION float
  linear(float r)
  {
    return r >= 1 || r <= -1? 0: 1 - Kokkos::Experimental::abs(r);
  }

  KOKKOS_INLINE_FUNCTION float
  derivative_linear(float r, float inv_cellsize)
  {
    return r >= 1 || r <= -1 || r == 0? 0: -Kokkos::Experimental::copysign(inv_cellsize, r);
  }

  KOKKOS_INLINE_FUNCTION float
  cubic_spline(float r, int ntype)
  {
    if (r >= 1 && r < 2)
      return ntype == 1? 0: ((-r/6 + 1)*r - 2)*r + 4.0/3;

    if (r >= 0 && r < 1)
      return ntype == -2 || ntype == 2? (r*r/6 - 1)*r + 1:
             ntype ==  1? (r/3 - 1)*r*r + 2.0/3:
                          (r/2 - 1)*r*r + 2.0/3;
    
    if (r >= -1 && r < 0)
      return ntype ==  -2 || ntype == 2? (-r*r/6 + 1)*r + 1:
             ntype == -1? (-r/3 - 1)*r*r + 2.0/3:
                          (-0.5*r - 1)*r*r + 2.0/3.0;
    
    if (r >= -2 && r < -1)
      return ntype == -1? 0 :((r/6 + 1)*r + 2)*r + 4.0/3;

    return 0;
  }

  KOKKOS_INLINE_FUNCTION float
  derivative_cubic_spline(float r, int ntype, float inv_cellsize)
  {
    if (r >= 1 && r < 2)
      return ntype == 1? 0 : inv_cellsize*((-0.5*r + 2)*r - 2);
    
    if (r >=0 && r < 1)
      return ntype == -2 || ntype == 2? inv_cellsize*(r*r/2 - 1):
             ntype ==  1? inv_cellsize*r*(r - 2):
                          inv_cellsize*(r*3/2 - 2)*r;
    
    if (r >= -1 && r < 0)
      return ntype == -2 || ntype == 2? inv_cellsize*(-r*r/2 + 1):
             ntype == -1? inv_cellsize*(-r - 2)*r:
                          inv_cellsize*(-r*3/2 - 2)*r;
    
    if (r >= -2 && r < -1)
      return ntype == -1? 0 : inv_cellsize*((r/2 + 2)*r + 2);

    return 0;
  }

  KOKKOS_INLINE_FUNCTION float
  bernstein_quadratic(float r, int ntype)
  {
    if (r >= 1.0 || r <= -1.0)
      return 0;

    if (ntype==1)
      return r > 0.5 || r < -0.5? 0: 0.5 - 2*r*r;

    r = Kokkos::Experimental::abs(r);
    return (1 - r)*(1-r);
  }

  KOKKOS_INLINE_FUNCTION float
  derivative_bernstein_quadratic(float r, int ntype, float inv_cellsize)
  {
    if (r >= 1 || r <= -1)
      return 0;

    if (ntype == 1)
      return r > 0.5 || r < -0.5? 0: -4*r*inv_cellsize;
    
	  return 2*(r - Kokkos::Experimental::copysign(1, r))*inv_cellsize;
  }

  KOKKOS_INLINE_FUNCTION float
  quadratic_spline(float r, int ntype)
  {
    switch (ntype)
    {
      case 0:
        return r >= -1.5 && r < -0.5? (r + 3)/2*r + 1.125:
               r >= -0.5 && r <  0.5? -r       *r + 0.75:
               r >=  0.5 && r <  1.5? (r - 3)/2*r + 1.125:
                                       0;
      case  2:
      case -2:
        return r >= -1.5 && r < -0.5? (0.5 * r + 1.5) * r + 1.125:
               r >= -0.5 && r < 0   ? 1 + r:
	       r >= 0    && r < 0.5 ? -r            + 1:
               r >= 0.5  && r < 1.5 ? (r/2 - 1.5)*r + 1.125:
                                    0;

      case -1:
        return r >= -1 && r < -0.5? 1 + r:
               r >= -0.5 && r < 0.5? -r * r + 0.75:
               r >= 0.5 && r < 1.5? (0.5 * r - 1.5) * r + 1.125:
                                    0;

      case 1:
        return r >= -1.5 && r < -0.5? (0.5 * r + 1.5) * r + 1.125:
               r >= -0.5 && r < 0.5? -r * r + 0.75:
               r >= 0.5 && r < 1? 1 - r:
                                  0;

      default:
        return 0;
    }
  }

  KOKKOS_INLINE_FUNCTION float
  derivative_quadratic_spline(float r, int ntype, float inv_cellsize)
  {
    switch (ntype)
    {
      case 0:
        return r >=  0.5 && r < 1.5? inv_cellsize*(r - 1.5):
               r >= -0.5 && r < 0.5? -2*inv_cellsize*r:
               r >= -1.5 && r < 0.5? inv_cellsize*(r + 1.5):
                                     0;
      case  2:
      case -2:
        return r >= -1.5 && r < -0.5? inv_cellsize*(r + 1.5):
               r >= -0.5 && r < 0   ? inv_cellsize:
	       r >= 0    && r < 0.5 ? -inv_cellsize:
               r >= 0.5  && r < 1.5 ? inv_cellsize*(r - 1.5):
                                    0;
      case -1:
        return r >= -1 && r < -0.5? inv_cellsize:
               r >= -0.5 && r < 0.5? -2 * inv_cellsize * r:
               r >= 0.5 && r < 1.5? inv_cellsize*(r - 1.5):
                                    0;
    
      case 1:
        return r >= -1.5 && r < -0.5? inv_cellsize*(r + 1.5):
               r >= -0.5 && r < 0.5? - 2 * inv_cellsize * r:
               r >= 0.5 && r < 1? -inv_cellsize:
                                  0;

      default:
        return 0;
    }
  }
}
#endif
