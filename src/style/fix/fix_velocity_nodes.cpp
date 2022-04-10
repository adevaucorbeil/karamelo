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
    xvalue = input->parsev(args[3]);

    string previous = args[3];

    // Replace "time" by "time - dt" in the x argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    xprevvalue = input->parsev(previous);

    v[0] = &input->expressions[args[3]];
    v_prev[0] = &input->expressions[previous];
  }
  else
    v_prev[0] = v[0] = nullptr;

  if (domain->dimension >= 2 && args[4] != "NULL")
  {
    yvalue = input->parsev(args[4]);
    yset = true;

    string previous = args[4];

    // Replace "time" by "time - dt" in the y argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    yprevvalue = input->parsev(previous);

    v[1] = &input->expressions[args[4]];
    v_prev[1] = &input->expressions[previous];
  }
  else
    v_prev[1] = v[1] = nullptr;

  if (domain->dimension == 3 && args[5] != "NULL")
  {
    zvalue = input->parsev(args[5]);
    zset = true;

    string previous = args[5];

    // Replace "time" by "time - dt" in the z argument:
    while(previous.find(time)!=std::string::npos)
      previous.replace(previous.find(time),time.length(),"time - dt");

    zprevvalue = input->parsev(previous);

    v[2] = &input->expressions[args[5]];
    v_prev[2] = &input->expressions[previous];
  }
  else
    v_prev[2] = v[2] = nullptr;
}

void FixVelocityNodes::prepare()
{
  xvalue.result(mpm);
  yvalue.result(mpm);
  zvalue.result(mpm);
  xprevvalue.result(mpm);
  yprevvalue.result(mpm);
  zprevvalue.result(mpm);

  ftot = Vector3d();
}

void FixVelocityNodes::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixVelocityNodes::post_update_grid_state(Grid &grid)
{
  // cout << "In FixVelocityNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:

  for (int i = 0; i < 3; i++)
    if (v[i])
      v[i]->evaluate(grid);

  for (int in = 0; in < grid.nnodes_local + grid.nnodes_ghost; in++)
  {
    if (!(grid.mask[in] & groupbit))
      continue;

	(*input->vars)["x"] = Var("x",   grid.x [in][0]);
	(*input->vars)["y"] = Var("y",   grid.x [in][1]);
	(*input->vars)["z"] = Var("z",   grid.x [in][2]);
	(*input->vars)["x0"] = Var("x0", grid.x0[in][0]);
	(*input->vars)["y0"] = Var("y0", grid.x0[in][1]);
	(*input->vars)["z0"] = Var("z0", grid.x0[in][2]);

    Vector3d Dv;
    if (xset)
    {
      double vx = xvalue.result(mpm, true);

      if (vx)
        cout << in << ": " << vx << " " << (*v[0])[in] << endl;

      double vx_old = xprevvalue.result(mpm, true);
      // cout << "Set v_update[0] to " << xvalue.eq() << "=" << vx << endl;
      // cout << "Set v[0] to " << vx_old << endl;
      Dv[0] = vx - grid.v_update[in][0];
      grid.v_update[in][0] = vx;
      grid.v[in][0] = vx_old;
    }
    if (yset)
    {
      double vy = yvalue.result(mpm, true);
      double vy_old = yprevvalue.result(mpm, true);
      // cout << "Set v_update[1] to " << "=" <<  vy << endl;
      // cout << "Set v[1] to " << "=" <<  vy_old << endl;
      Dv[1] = vy - grid.v_update[in][1];
      grid.v_update[in][1] = vy;
      grid.v[in][1] = vy_old;
    }
    if (zset)
    {
      double vz = zvalue.result(mpm, true);
      double vz_old = zprevvalue.result(mpm, true);
      // cout << "Set v_update[2] to " << "=" <<  vz << endl;
      // cout << "Set v[2] to " << "=" <<  vz_old << endl;
      Dv[2] = vz - grid.v_update[in][2];
      grid.v_update[in][2] = vz;
      grid.v[in][2] = vz_old;
    }
    ftot += grid.mass[in]*Dv/update->dt;
  }
}

void FixVelocityNodes::post_velocities_to_grid(Grid &grid) {
  // cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  for (int in = 0; in < grid.nnodes_local + grid.nnodes_ghost; in++)
  {
    if (!(grid.mask[in] & groupbit))
      continue;

    if (xset) grid.v[in][0] = xvalue.result(mpm, true);
    if (yset) grid.v[in][1] = yvalue.result(mpm, true);
    if (zset) grid.v[in][2] = zvalue.result(mpm, true);
  }
}

void FixVelocityNodes::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.write_to_restart(of);
    xprevvalue.write_to_restart(of);
  }
  if (yset) {
    yvalue.write_to_restart(of);
    yprevvalue.write_to_restart(of);
  }
  if (zset) {
    zvalue.write_to_restart(of);
    zprevvalue.write_to_restart(of);
  }
}

void FixVelocityNodes::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
   ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.read_from_restart(ifr);
    xprevvalue.read_from_restart(ifr);
  }
  if (yset) {
    yvalue.read_from_restart(ifr);
    yprevvalue.read_from_restart(ifr);
  }
  if (zset) {
    zvalue.read_from_restart(ifr);
    zprevvalue.read_from_restart(ifr);
  }
}
