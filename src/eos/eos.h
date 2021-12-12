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

#ifndef MPM_EOS_H
#define MPM_EOS_H

#include <pointers.h>
#include <vector>
#include <matrix.h>

/*! Parent class of all the different kinds of EOS (equations of state) that can be used.
 *
 * Stores the EOS id, the functions that returns the bulk modulus and the reference density
 * and the method to compute the pressure.
 */
class EOS : protected Pointers {
 public:
  string id;                        ///< EOS identification string
  string style;                     ///< EOS style

  EOS(class MPM *, vector<string>);
  virtual ~EOS();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each EOS
  //virtual compute_pressure()
  virtual double rho0() = 0;
  virtual double K() = 0;
  virtual void compute_pressure(double &, double &, const double, const double, const double, const Matrix3d, const double, const double T = 0) = 0;

  virtual void write_restart(ofstream*) = 0;
  virtual void read_restart(ifstream*) = 0;
  //protected:
};

#endif

/*! \defgroup eos eos

\section Syntax Syntax
\code
eos(eos-ID, eos_type, eos_specific_arguments)
\endcode

<ul>
<li>eos-ID: name of the eos to be created.</li>
<li>eos-type: linear, shock, fluid ...</li>
<li>eos_specific_arguments: list of arguments specific to the eos type used. </li>
</ul>

\section Examples Examples
\code
eos(eosl, linear, rho, K)
\endcode
Defines a linear equation of state (EOSLinear) called 'eosl' with a reference bulk density rho, a bulk modulus K

\code
eos(eoss, shock, rho0, K, c0, S, Gamma, cv, Tr, Q1, Q2)
\endcode
Defines a Mie-Gruneisen equation of state (EOSShock) called 'eoss' with a reference bulk density rho, a bulk modulus K, a bulk speed of sound c0,
a Hugoniot slope coefficient S, a Grüneisen Gamma at the reference state Gamma, a heat capacity at constant
volume cv, temperature of reference Tr, and artificial viscosity coefficients Q1 and Q2.


\section Description Description

This command defines an equation of state used to determine the hydrostatic pressure.

\section EOS_Subclasses EOS types

To access the different eos type supported, and how to use the eos() command with them, refer to the corresponding class.
*/

