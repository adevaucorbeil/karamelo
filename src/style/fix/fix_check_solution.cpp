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

using namespace std;
using namespace FixConst;
using namespace MathSpecial;


FixChecksolution::FixChecksolution(MPM *mpm, vector<string> args):
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

  if (domain->dimension == 3 && args.size()<6) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 2 && args.size()<5) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 5 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR,"Error: too few arguments for fix_check_solution: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"_check_solution needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixChecksolution with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (domain->dimension >= 2) {
    if (args[4] != "NULL") {
      yvalue = input->parsev(args[4]);
      yset = true;
    }
  }

  if (domain->dimension == 3) {
    if (args[5] != "NULL") {
      zvalue = input->parsev(args[5]);
      zset = true;
    }
  }
}

void FixChecksolution::prepare()
{
  xvalue.result(mpm);
  yvalue.result(mpm);
  zvalue.result(mpm);

  error_vec = Vector3d();
  u_th = Vector3d();
}

void FixChecksolution::reduce()
{
  Vector3d error_reduced, u_th_reduced;

  // Reduce error:
  MPI_Allreduce(error_vec.elements, error_reduced.elements, 3, MPI_DOUBLE,MPI_SUM, universe->uworld);
  MPI_Allreduce(u_th     .elements, u_th_reduced .elements, 3, MPI_DOUBLE,MPI_SUM, universe->uworld);

  double vtot = 0;

  int solid = group->solid[igroup];

  if (solid == -1)
    for (const Solid *solid: domain->solids)
      vtot += solid->vtot;
  else
    vtot += domain->solids.at(solid)->vtot;

  (*input->vars)[id + "_s"] = Var(id + "_s", sqrt((error_reduced[0] + error_reduced[1] + error_reduced[2])/vtot));
  (*input->vars)[id + "_x"] = Var(id + "_x", (*input->vars)[id + "_x"].result() + update->dt*(error_reduced[0] + error_reduced[1] + error_reduced[2]));
  (*input->vars)[id + "_y"] = Var(id + "_y", (*input->vars)[id + "_y"].result() + update->dt*(u_th_reduced[0] + u_th_reduced[1] + u_th_reduced[2]));
  (*input->vars)[id + "_z"] = Var(id + "_z", sqrt((*input->vars)[id + "_x"].result()/(*input->vars)[id + "_y"].result()));
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}

void FixChecksolution::final_integrate(Solid &solid, int ip)
{
  if ((update->ntimestep != output->next && update->ntimestep != update->nsteps) ||
      !(solid.mask.at(ip) & groupbit))
    return;

  (*input->vars)["x0"] = Var("x0", solid.x0.at(ip)[0]);
  (*input->vars)["y0"] = Var("y0", solid.x0.at(ip)[1]);
  (*input->vars)["z0"] = Var("z0", solid.x0.at(ip)[2]);

  if (xset)
  {
    double ux = xvalue.result(mpm, true);
    error_vec[0] += solid.vol0.at(ip)*square(ux - (solid.x[ip][0] - solid.x0.at(ip)[0]));
    u_th[0] += solid.vol0.at(ip)*ux*ux;                  
  }                                                
  if (yset)
  {                                      
    double uy = yvalue.result(mpm, true);                       
    error_vec[1] += solid.vol0.at(ip)*square(uy - (solid.x[ip][1] - solid.x0.at(ip)[1]));
    u_th[1] += solid.vol0.at(ip)*uy*uy;                  
  }                                                
  if (zset)
  {                                      
    double uz = zvalue.result(mpm, true);                       
    error_vec[2] += solid.vol0.at(ip)*square(uz - (solid.x[ip][2] - solid.x0.at(ip)[2]));
    u_th[2] += solid.vol0.at(ip)*uz*uz;
  }
}


void FixChecksolution::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset) {
    string eq = xvalue.eq();
    size_t N = eq.size();
    double value = xvalue.result();
    bool cst = xvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }

  if (yset) {
    string eq = yvalue.eq();
    size_t N = eq.size();
    double value = yvalue.result();
    bool cst = yvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }

  if (zset) {
    string eq = zvalue.eq();
    size_t N = eq.size();
    double value = zvalue.result();
    bool cst = zvalue.is_constant();
    of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
    of->write(reinterpret_cast<const char *>(eq.c_str()), N);
    of->write(reinterpret_cast<const char *>(&value), sizeof(double));
    of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
  }
}

void FixChecksolution::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    xvalue = Var(eq, value, cst);
  }

  if (yset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    yvalue = Var(eq, value, cst);
  }

  if (zset) {
    string eq = "";
    size_t N = 0;
    double value = 0;
    bool cst = false;

    ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
    eq.resize(N);

    ifr->read(reinterpret_cast<char *>(&eq[0]), N);
    ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
    ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
    zvalue = Var(eq, value, cst);
  }
}
