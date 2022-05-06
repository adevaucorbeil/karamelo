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

#ifdef FIX_CLASS

FixStyle(check_solution,FixCheckSolution)

#else

#ifndef MPM_FIX_CHECK_SOLUTION_H
#define MPM_FIX_CHECK_SOLUTION_H

#include <fix.h>
#include <matrix.h>

class Expression;

/*! \ingroup fix fixchecksolution fix_check_solution

\section Syntax Syntax
\code
fix(fix-ID, check_solution, group-ID, ux, uy, uz)
\endcode


<ul>
<li>fix-ID: name of the fix to be created.</li>
<li>fix-type: velocity_nodes, force_nodes, kinetic_energy, ...</li>
<li>group-ID: name of the group onto which the fix will be applied. 
If 'all' is used, all particles will be selected.</li>
<li>ux, uy and uz: solution of the particle's displacements in the x, y and z directions. </li>
</ul>

\section Examples Examples
\code
c = sqrt(E/rho)

G = 0.05
phi2 = PI
phi3 = PI

ux = G*sin(PI*x)*sin(c*PI*time)
uy = G*sin(PI*y)*sin(c*PI*time + phi2)
uz = 0
fix(error, check_solution, all, ux, uy, uz)
\endcode
This fix calculates the error between all the particles' current position and their theoretical position given by \f$\boldsymbol{u}=[u_x(x,t), u_y(y,t), u_z]\f$. The resulting error can be accessed using the error_s, error_x, error_y and error_z variables.

\section Description Description

This command determine the error in displacement of the particles from the group 'group-ID', being provided the solution \f$\boldsymbol{u}=[ux, uy, uz]\f$, where ux, uy and uz are functions of positions and time. The fix generates four variables:
<ul>
<li>fix-ID_s:
\f{equation}{
  e(t^n) = \sqrt{\frac{\sum_{p=0}^{N_p} V_p^0\left|\left|\boldsymbol{u}_p^h(t^n) - \boldsymbol{u}(\boldsymbol{X}_p, t^n)\right|\right|^2}{V_{tot}^0}}
\f}</li>
<li>fix-ID_x:
\f{equation}{
  e = \sum_{t^n=0}^{T} \sum_{p=0}^{N_p} V_p^0\left|\left|\boldsymbol{u}_p^h(t^n) - \boldsymbol{u}(\boldsymbol{X}_p, t^n)\right|\right|^2
\f}</li>
<li>fix-ID_y:
\f{equation}{
  e = \sum_{t^n=0}^{T} \sum_{p=0}^{N_p} V_p^0\left|\left|\boldsymbol{u}(\boldsymbol{X}_p, t^n)\right|\right|^2
\f}</li>
<li>fix-ID_z:
\f{equation}{
  e_2 = \sqrt{\frac{\sum_{t^n=0}^{T} \sum_{p=0}^{N_p} V_p^0\left|\left|\boldsymbol{u}_p^h(t^n) - \boldsymbol{u}(\boldsymbol{X}_p, t^n)\right|\right|^2}{\sum_{t^n=0}^{T} \sum_{p=0}^{N_p} V_p^0\left|\left|\boldsymbol{u}(\boldsymbol{X}_p, t^n)\right|\right|^2}}
\f}</li>

\section Class Class description
*/
class FixCheckSolution : public Fix {
 public:
  FixCheckSolution(MPM *, vector<string>);

  void prepare();
  void reduce();

  void final_integrate(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, check_solution, group, ux)\n"},
      {2, "Usage: fix(fix-ID, check_solution, group, ux, uy)\n"},
      {3, "Usage: fix(fix-ID, check_solution, group, ux, uy, uz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  bool xset, yset, zset;               // Does the fix set the x, y, and z forces of the group?
  Vector3d error_vec, u_th;
  Expression *u[3];
};

#endif
#endif

