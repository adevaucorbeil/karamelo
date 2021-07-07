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


using namespace std;

namespace BasisFunction { 
  inline double linear(const double r_, const int ntype)
  {
    double r = fabs(r_);
    if (r >= 1.0)
      return 0.0;
    else
      return 1.0 - r;
  }

  inline double derivative_linear(const double r, const int ntype, const double inv_cellsize)
  {
    if (r >= 1.0 || r <= -1.0 || r == 0)
      return 0.0;
    else if (r > 0.0)
      return -inv_cellsize;
    else
      return inv_cellsize;
  }

  inline double cubic_spline(const double r, const int ntype)
  {
    if (r >= 1 && r < 2) {
      if (ntype==1) {
	return 0;
      } else {
	return ((-1.0/6.0*r + 1)*r - 2)*r + 4.0/3.0;
      }
    } else if (r >=0 && r < 1) {
      if (ntype==-2) {
	return  (1.0/6.0*r*r - 1)*r + 1;
      } else if (ntype==2) {
	return 1;
      } else if (ntype==1) {
	return (1.0/3.0*r - 1)*r*r + 2.0/3.0;
      } else {
	return (0.5*r - 1)*r*r + 2.0/3.0;
      }
    } else if (r >= -1 && r < 0) {
      if (ntype==2) {
	return (-1.0/6.0*r*r + 1)*r + 1;
      } else if (ntype==-1) {
	return  (-1.0/3.0*r - 1)*r*r + 2.0/3.0;
      } else {
	return (-0.5*r - 1)*r*r + 2.0/3.0;
      }
    } else if (r >= -2 && r < -1) {
      return ((1.0/6.0*r + 1)*r + 2)*r + 4.0/3.0;
    } else {
      return 0;
    }
  }

  inline double derivative_cubic_spline(const double r, const int ntype, const double icellsize)
  {
    if (r >= 1 && r < 2) {
      if (ntype==1) {
	return -icellsize;// * (-1);
      } else {
	return icellsize * ((-0.5*r + 2)*r - 2);
      }
    } else if (r >=0 && r < 1) {
      if (ntype==-2) {
        return icellsize * (0.5 * r * r - 1);
      } else if (ntype==2) {
	return icellsize;// * (1);
      } else if (ntype==1) {
        return icellsize * r * (r - 2);
      } else {
        return icellsize * (3.0 / 2.0 * r - 2) * r;
      }
    } else if (r >= -1 && r < 0) {
      if (ntype==2) {
        return icellsize * (-0.5 * r * r + 1);
      } else if (ntype==-1) {
	return icellsize * (-r - 2)*r;
      } else {
	return icellsize * (-3.0/2.0*r - 2)*r;
      }
    } else if (r >= -2 && r < -1) {
      return icellsize * ((0.5 * r + 2) * r + 2);
    } else {
      return 0;
    }
  }

  inline double bernstein_quadratic(const double r_, const int ntype)
  {
    double r = fabs(r_);
    if (r >= 1.0) return 0;

    if (ntype==1) {
      // Inside node:
      if (r >= 0.5) return 0;
      else return 0.5-2*r*r;
    } else {
      // Edge node:
      return (1-r)*(1-r);
    }
  }

  inline double derivative_bernstein_quadratic(const double r_signed, const int ntype, const double icellsize)
  {
    double r = fabs(r_signed);
    if (r >= 1.0) return 0;

    if (ntype==1) {
      // Inside node:
      if (r > 0.5) return 0;
      return -4*r_signed*icellsize;
    } else {
      // Edge node:
      if (r_signed>0) {
	return -2*(1-r_signed)*icellsize;
      } else {
	return 2*(1+r_signed)*icellsize;
      }
    }
  }

  inline double quadratic_spline(const double r, const int ntype) {
    if (ntype == 0) {
      if (r >= 0.5 && r < 1.5) {
        return (0.5 * r - 1.5) * r + 1.125;
      } else if (r >= -0.5 && r < 0.5) {
        return -r * r + 0.75;
      } else if (r >= -1.5 && r < 0.5) {
        return (0.5 * r + 1.5) * r + 1.125;
      } else {
        return 0;
      }
    } else if (ntype == -2) {
      if (r >= 0. && r < 0.5) {
        return 1 - r;
      } else if (r >= 0.5 && r < 1.5) {
        return (0.5 * r - 1.5) * r + 1.125;
      } else {
        return 0;
      }
    } else if (ntype == -1) {
      if (r >= -1. && r < -0.5) {
        return 1 + r;
      } else if (r >= -0.5 && r < 0.5) {
        return -r * r + 0.75;
      } else if (r >= 0.5 && r < 1.5) {
        return (0.5 * r - 1.5) * r + 1.125;
      } else {
        return 0;
      }
    } else if (ntype == 1) {
      if (r >= -1.5 && r < -0.5) {
        return (0.5 * r + 1.5) * r + 1.125;
      } else if (r >= -0.5 && r < 0.5) {
        return -r * r + 0.75;
      } else if (r >= 0.5 && r < 1.) {
        return 1 - r;
      } else {
        return 0;
      }
    } else {
      //ntype == 2
      if (r >= -1.5 && r < -0.5) {
        return (0.5 * r + 1.5) * r + 1.125;
      } else if (r >= -0.5 && r <= 0.) {
        return 1 + r;
      } else {
        return 0;
      }
    }
  }

  inline double derivative_quadratic_spline(const double r, const int ntype,
                                            const double icellsize) {
    if (ntype == 0) {
      if (r >= 0.5 && r < 1.5) {
        return icellsize*(r - 1.5);
      } else if (r >= -0.5 && r < 0.5) {
        return - 2 * icellsize * r;
      } else if (r >= -1.5 && r < 0.5) {
        return icellsize*(r + 1.5);
      } else {
        return 0;
      }
    } else if (ntype == -2) {
      if (r >= 0. && r < 0.5) {
        return -icellsize;
      } else if (r >= 0.5 && r < 1.5) {
        return icellsize*(r - 1.5);
      } else {
        return 0;
      }
    } else if (ntype == -1) {
      if (r >= -1. && r < -0.5) {
        return icellsize;
      } else if (r >= -0.5 && r < 0.5) {
        return -2 * icellsize * r;
      } else if (r >= 0.5 && r < 1.5) {
        return icellsize*(r - 1.5);
      } else {
        return 0;
      }
    } else if (ntype == 1) {
      if (r >= -1.5 && r < -0.5) {
        return icellsize*(r + 1.5);
      } else if (r >= -0.5 && r < 0.5) {
        return - 2 * icellsize * r;
      } else if (r >= 0.5 && r < 1.) {
        return -icellsize;
      } else {
        return 0;
      }
    } else {
      // ntype == 2
      if (r >= -1.5 && r < -0.5) {
        return icellsize*(r + 1.5);
      } else if (r >= -0.5 && r <= 0.) {
        return icellsize;
      } else {
        return 0;
      }
    }
  }
}
#endif
