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


#include <fix_velocity_nodes.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <grid.h>
#include <error.h>
#include <update.h>
#include <universe.h>
#include <expression_operation.h>

using namespace std;
using namespace FixConst;


FixVelocityNodes::FixVelocityNodes(MPM *mpm, vector<string> args):
  Fix(mpm, args, POST_UPDATE_GRID_STATE | POST_VELOCITIES_TO_GRID)
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

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->one(FLERR, "fix_velocity_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of "+ group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixVelocityNodes with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  string time = "time";

  if (args[3] != "NULL")
  {
    // xvalue = input->parsev(args[3]);
    xset = true;
    input->parsev(args[3]);

    string previous = args[3];

    // Replace "time" by "time - dt" in the x argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    input->parsev(previous);

    v[0] = &input->expressions[args[3]];
    v_prev[0] = &input->expressions[previous];
  }
  else
    v_prev[0] = v[0] = nullptr;

  if (domain->dimension >= 2 && args[4] != "NULL")
  {
    input->parsev(args[4]);
    yset = true;

    string previous = args[4];

    // Replace "time" by "time - dt" in the y argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    input->parsev(previous);

    v[1] = &input->expressions[args[4]];
    v_prev[1] = &input->expressions[previous];
  }
  else
    v_prev[1] = v[1] = nullptr;

  if (domain->dimension == 3 && args[5] != "NULL")
  {
    input->parsev(args[5]);
    zset = true;

    string previous = args[5];

    // Replace "time" by "time - dt" in the z argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    input->parsev(previous);

    v[2] = &input->expressions[args[5]];
    v_prev[2] = &input->expressions[previous];
  }
  else
    v_prev[2] = v[2] = nullptr;
}

void FixVelocityNodes::prepare()
{
  ftot = Vector3d();
}

void FixVelocityNodes::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  //(*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  //(*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  //(*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixVelocityNodes::post_update_grid_state(Grid &grid)
{
  // cout << "In FixVelocityNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:

  for (int i = 0; i < 3; i++)
    if (v[i])
      v[i]->evaluate(grid);

  int groupbit = this->groupbit;
  double dt = update->dt;
  Kokkos::View<double*> mass = grid.mass;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<Vector3d*> gv = grid.v, v_update = grid.v_update;

  for (int i = 0; i < 3; i++)
    if (v[i])
    {
      Kokkos::View<double **> v_i = v[i]->registers;
      Kokkos::View<double **> v_prev_i = v_prev[i]->registers;

      Kokkos::parallel_reduce("FixVelocityNodes::post_update_grid_state", grid.nnodes_local + grid.nnodes_ghost,
                              KOKKOS_LAMBDA(const int &in, double &ftot_i)
      {
        if (!(mask[in] & groupbit))
          return;

        gv[in][i] = v_prev_i(0, in);

        ftot_i += mass[in]*((v_update[in][i] = v_i(0, in)) - v_update[in][i])/dt;
      }, ftot[i]);
    }
}

void FixVelocityNodes::post_velocities_to_grid(Grid &grid) {
  // cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  /*for (int in = 0; in < grid.nnodes_local + grid.nnodes_ghost; in++)
  {
    if (!(grid.mask[in] & groupbit))
      continue;

    if (xset) grid.v[in][0] = (*v[0])[in];
    if (yset) grid.v[in][1] = (*v[1])[in];
    if (zset) grid.v[in][2] = (*v[2])[in];
  }*/
  error->all(FLERR, "FixVelocityNodes::post_velocities_to_grid not parallelized");
}

void FixVelocityNodes::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  //if (xset) {
  //  xvalue.write_to_restart(of);
  //  xprevvalue.write_to_restart(of);
  //}
  //if (yset) {
  //  yvalue.write_to_restart(of);
  //  yprevvalue.write_to_restart(of);
  //}
  //if (zset) {
  //  zvalue.write_to_restart(of);
  //  zprevvalue.write_to_restart(of);
  //}
}

void FixVelocityNodes::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
   ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  //if (xset) {
  //  xvalue.read_from_restart(ifr);
  //  xprevvalue.read_from_restart(ifr);
  //}
  //if (yset) {
  //  yvalue.read_from_restart(ifr);
  //  yprevvalue.read_from_restart(ifr);
  //}
  //if (zset) {
  //  zvalue.read_from_restart(ifr);
  //  zprevvalue.read_from_restart(ifr);
  //}
}
