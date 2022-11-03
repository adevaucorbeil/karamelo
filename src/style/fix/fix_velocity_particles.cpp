/* ----------------------------------------------------------------------
*
*                    ***       Karamelo       ***
*               Parallel Material Point Method Simulator
* 
* Copyright (2020) Alban de Vaucorbeil, alban.devaucorbeil@deakin.edu.au
* Institute for Frontier Materials, Deakin University
* Geelong VIC 3216, Australia

* This software is distributed under the GNU General Public License.
*
* ----------------------------------------------------------------------- */

#include <fix_velocity_particles.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <grid.h>
#include <error.h>
#include <update.h>
#include <method.h>
#include <universe.h>
#include <expression_operation.h>


using namespace std;
using namespace FixConst;



FixVelocityParticles::FixVelocityParticles(MPM *mpm, vector<string> args):
 Fix(mpm, args, INITIAL_INTEGRATE | POST_ADVANCE_PARTICLES)
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
   return;
 }

 if (args.size() < Nargs.find(domain->dimension)->second) {
   error->all(FLERR, "Error: too few arguments for fix_velocity_nodes.\n" +
                         usage.find(domain->dimension)->second);
 }

 if (group->pon[igroup].compare("particles") !=0 ) {
   error->one(FLERR, "fix_velocity_nodes needs to be given a group of nodes" +
                         group->pon[igroup] + ", " + args[2] +
                         " is a group of " + group->pon[igroup] + ".\n");
 }
 if (universe->me == 0) {
   cout << "Creating new fix FixVelocityParticles with ID: " << args[0] << endl;
 }
 id = args[0];

 v_prev[0] = v[0] = nullptr;
 v_prev[1] = v[1] = nullptr;
 v_prev[2] = v[2] = nullptr;


 string time = "time";

  for (int dim = 0; dim < domain->dimension; dim++) {
    if (args[3 + dim] != "NULL") {

      input->parsev(args[3 + dim]);

      string previous = args[3 + dim];
      // Replace "time" by "time - dt" in the x argument:
      while (previous.find(time) != std::string::npos)
        previous.replace(previous.find(time), time.length(), "time - dt");

      input->parsev(previous);

      v[dim] = &input->expressions[args[3 + dim]];
      v_prev[dim] = &input->expressions[previous];
    }
  }
}

void FixVelocityParticles::prepare()
{
 ftot = Vector3d();
}

void FixVelocityParticles::reduce()
{
 Vector3d ftot_reduced;

 // Reduce ftot:
 MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM,
               universe->uworld);

 input->parsev(id + "_x", ftot_reduced[0]);
 input->parsev(id + "_y", ftot_reduced[1]);
 input->parsev(id + "_z", ftot_reduced[2]);

 // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
 // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
 // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixVelocityParticles::initial_integrate(Solid &solid)
{
 // Go through all the particles in the group and set v_update to the right value:

  for (int i = 0; i < 3; i++) {
    if (v[i])
      v[i]->evaluate(solid);
    if (v_prev[i])
      v_prev[i]->evaluate(solid);
  }

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<Vector3d*> sv = solid.v, sv_update = solid.v_update;

  for (int i = 0; i < 3; i++)
    if (v[i])
    {
      Kokkos::View<float **> v_i = v[i]->registers;
      Kokkos::View<float **> v_prev_i = v_prev[i]->registers;


      Kokkos::parallel_for("FixVelocityParticles::initial_integrate", solid.np_local,
      KOKKOS_LAMBDA(const int &ip)
      {
        if (!(mask[ip] & groupbit))
          return;

	sv_update[ip][i] = v_i(0, ip);
	sv[ip][i] = v_prev_i(0, ip);
      });
    }
}

void FixVelocityParticles::post_advance_particles(Solid &solid) {
 // Go through all the particles in the group and set v to the right value:

  for (int i = 0; i < 3; i++) {
    if (v[i])
      v[i]->evaluate(solid);
  }

  int groupbit = this->groupbit;
  float dt = update->dt;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<float*> smass = solid.mass;
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<Vector3d*> sv = solid.v, sv_update = solid.v_update;

  for (int i = 0; i < 3; i++)
    if (v[i])
    {
      Kokkos::View<float **> v_i = v[i]->registers;

      Kokkos::parallel_reduce("FixVelocityParticles::post_advance_particles", solid.np_local,
      KOKKOS_LAMBDA(const int &ip, float &ftot_i)
      {
        if (!(mask[ip] & groupbit))
          return;

	const float &xold_i = sx[ip][i] - dt*sv_update[ip][i];
	const float &Dv_i   = v_i(0, ip) - sv[ip][i];

	sv[ip][i] = v_i(0, ip);
	sx[ip][i] = xold_i + dt * v_i(0, ip);
	ftot_i += smass[ip] * Dv_i / dt;
      }, ftot[i]);
    }
}

void FixVelocityParticles::write_restart(ofstream *of) {
 // of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
 // of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
 // of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

 // if (xset) {
 //   xvalue.write_to_restart(of);
 //   xprevvalue.write_to_restart(of);
 // }
 // if (yset) {
 //   yvalue.write_to_restart(of);
 //   yprevvalue.write_to_restart(of);
 // }
 // if (zset) {
 //   zvalue.write_to_restart(of);
 //   zprevvalue.write_to_restart(of);
 // }
}

void FixVelocityParticles::read_restart(ifstream *ifr) {
//  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
//  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
//  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

//  if (xset) {
//    xvalue.read_from_restart(ifr);
//    xprevvalue.read_from_restart(ifr);
//  }
//  if (yset) {
//    yvalue.read_from_restart(ifr);
//    yprevvalue.read_from_restart(ifr);
//  }
//  if (zset) {
//    zvalue.read_from_restart(ifr);
//    zprevvalue.read_from_restart(ifr);
//  }
}
