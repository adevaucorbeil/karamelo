/* ----------------------------------------------------------------------
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


#include <fix_check_solution.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <input.h>
#include <update.h>
#include <output.h>
#include <math_special.h>
#include <universe.h>
#include <solid.h>
#include <error.h>
#include <var.h>
#include <expression_operation.h>

using namespace std;
using namespace FixConst;
using namespace MathSpecial;


FixCheckSolution::FixCheckSolution(MPM *mpm, vector<string> args):
  Fix(mpm, args, FINAL_INTEGRATE)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];

    xset = yset = zset = false;
    return;
  }

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_velocity_nodes.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"_check_solution needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixCheckSolution with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    input->parsev(args[3]);
    xset = true;
    u[0] = &input->expressions[args[3]];
  }
  else
    u[0] = nullptr;

  if (domain->dimension >= 2 && args[4] != "NULL") {
    input->parsev(args[4]);
    yset = true;
    u[1] = &input->expressions[args[4]];
  }
  else
    u[1] = nullptr;

  if (domain->dimension == 3 && args[5] != "NULL") {
    input->parsev(args[5]);
    zset = true;
    u[2] = &input->expressions[args[5]];    
  } else
    u[2] = nullptr;
}

void FixCheckSolution::prepare()
{
  error_vec = Vector3d();
  u_th = Vector3d();
}

void FixCheckSolution::reduce()
{
  Vector3d error_reduced, u_th_reduced;

  // Reduce error:
  MPI_Allreduce(error_vec.elements, error_reduced.elements, 3, MPI_FLOAT,MPI_SUM, universe->uworld);
  MPI_Allreduce(u_th     .elements, u_th_reduced .elements, 3, MPI_FLOAT,MPI_SUM, universe->uworld);

  float vtot = 0;

  int solid = group->solid[igroup];

  if (solid == -1)
    for (const Solid *solid: domain->solids)
      vtot += solid->vtot;
  else
    vtot += domain->solids[solid]->vtot;

  input->parsev(id + "_s", sqrt((error_reduced[0] + error_reduced[1] + error_reduced[2])/vtot));
  input->parsev(id + "_x", (*input->vars)[id + "_x"].result() + update->dt*(error_reduced[0] + error_reduced[1] + error_reduced[2]));
  input->parsev(id + "_y", (*input->vars)[id + "_y"].result() + update->dt*(u_th_reduced[0] + u_th_reduced[1] + u_th_reduced[2]));
  input->parsev(id + "_z", sqrt((*input->vars)[id + "_x"].result()/(*input->vars)[id + "_y"].result()));
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}

void FixCheckSolution::final_integrate(Solid &solid)
{
  // Go through all the nodes in the group and set v to the right value:

  for (int i = 0; i < 3; i++)
    if (u[i])
      u[i]->evaluate(solid);

  int groupbit = this->groupbit;
  float ntimestep = update->ntimestep;
  int nsteps = update->nsteps;
  bigint next = output->next;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<float*> vol0 = solid.vol0;
  Kokkos::View<Vector3d*> x = solid.x;
  Kokkos::View<Vector3d*> x0 = solid.x0;

  for (int i = 0; i < 3; i++)
    if (u[i])
    {
      Kokkos::View<float **> u_i = u[i]->registers;

      Kokkos::parallel_reduce(
          "FixCheckSolution::final_integrate", solid.np_local,
          KOKKOS_LAMBDA(const int &ip, float &error_vec_i, float &u_th_i) {
            if ((ntimestep != next && ntimestep != nsteps) ||
                !(mask[ip] & groupbit))
              return;

            error_vec_i += vol0[ip] *
	      Kokkos::Experimental::pow(u_i(0, ip) - (x[ip][i] - x0[ip][i]), 2);

            u_th_i += vol0[ip] * u_i(0, ip) * u_i(0, ip);
          }, error_vec[i], u_th[i]);
    }
}


void FixCheckSolution::write_restart(ofstream *of) {
  error->all(FLERR, "Error: Restart broken.\n");
  // of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  // of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  // of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  // if (xset) {
  //   string eq = xvalue.eq();
  //   size_t N = eq.size();
  //   float value = xvalue.result();
  //   bool cst = xvalue.is_constant();
  //   of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  //   of->write(reinterpret_cast<const char *>(eq.c_str()), N);
  //   of->write(reinterpret_cast<const char *>(&value), sizeof(float));
  //   of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  // }

  // if (yset) {
  //   string eq = yvalue.eq();
  //   size_t N = eq.size();
  //   float value = yvalue.result();
  //   bool cst = yvalue.is_constant();
  //   of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  //   of->write(reinterpret_cast<const char *>(eq.c_str()), N);
  //   of->write(reinterpret_cast<const char *>(&value), sizeof(float));
  //   of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  // }

  // if (zset) {
  //   string eq = zvalue.eq();
  //   size_t N = eq.size();
  //   float value = zvalue.result();
  //   bool cst = zvalue.is_constant();
  //   of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  //   of->write(reinterpret_cast<const char *>(eq.c_str()), N);
  //   of->write(reinterpret_cast<const char *>(&value), sizeof(float));
  //   of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  // }
}

void FixCheckSolution::read_restart(ifstream *ifr) {
  error->all(FLERR, "Error: Restart broken.\n");
  // ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  // ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  // ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  // if (xset) {
  //   string eq = "";
  //   size_t N = 0;
  //   float value = 0;
  //   bool cst = false;

  //   ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  //   eq.resize(N);

  //   ifr->read(reinterpret_cast<char *>(&eq[0]), N);
  //   ifr->read(reinterpret_cast<char *>(&value), sizeof(float));
  //   ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
  //   xvalue = Var(eq, value, cst);
  // }

  // if (yset) {
  //   string eq = "";
  //   size_t N = 0;
  //   float value = 0;
  //   bool cst = false;

  //   ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  //   eq.resize(N);

  //   ifr->read(reinterpret_cast<char *>(&eq[0]), N);
  //   ifr->read(reinterpret_cast<char *>(&value), sizeof(float));
  //   ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
  //   yvalue = Var(eq, value, cst);
  // }

  // if (zset) {
  //   string eq = "";
  //   size_t N = 0;
  //   float value = 0;
  //   bool cst = false;

  //   ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  //   eq.resize(N);

  //   ifr->read(reinterpret_cast<char *>(&eq[0]), N);
  //   ifr->read(reinterpret_cast<char *>(&value), sizeof(float));
  //   ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
  //   zvalue = Var(eq, value, cst);
  // }
}
