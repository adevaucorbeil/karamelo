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

#include <fix_cutting_tool.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


// void write_var(Var v, ofstream *of) {
//   string eq = v.eq();
//   size_t N = eq.size();
//   double value = v.result();
//   bool cst = v.is_constant();
//   of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
//   of->write(reinterpret_cast<const char *>(eq.c_str()), N);
//   of->write(reinterpret_cast<const char *>(&value), sizeof(double));
//   of->write(reinterpret_cast<const char *>(&cst), sizeof(bool));
// }

// Var read_var(ifstream *ifr) {
//   string eq = "";
//   size_t N = 0;
//   double value = 0;
//   bool cst = false;

//   ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
//   eq.resize(N);

//   ifr->read(reinterpret_cast<char *>(&eq[0]), N);
//   ifr->read(reinterpret_cast<char *>(&value), sizeof(double));
//   ifr->read(reinterpret_cast<char *>(&cst), sizeof(bool));
//   return Var(eq, value, cst);
// }

FixCuttingTool::FixCuttingTool(MPM *mpm, vector<string> args)
    : Fix(mpm, args, INITIAL_INTEGRATE) {
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

    K = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR, "fix_cuttingtool needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixCuttingTool with ID: " << args[0] << endl;
  }
  id = args[0];
  K = input->parsev(args[3]).result(mpm);

  xtvalue = input->parsev(args[4]);
  ytvalue = input->parsev(args[5]);
  ztvalue = input->parsev(args[6]);

  vtxvalue = input->parsev(args[7]);
  vtyvalue = input->parsev(args[8]);
  vtzvalue = input->parsev(args[9]);

  xAvalue = input->parsev(args[10]);
  yAvalue = input->parsev(args[11]);

  xBvalue = input->parsev(args[12]);
  yBvalue = input->parsev(args[13]);
}

void FixCuttingTool::prepare()
{
  xtvalue .result(mpm);
  ytvalue .result(mpm);
  ztvalue .result(mpm);
  vtxvalue.result(mpm);
  vtyvalue.result(mpm);
  vtzvalue.result(mpm);
  xAvalue .result(mpm);
  yAvalue .result(mpm);
  xBvalue .result(mpm);
  yBvalue .result(mpm);

  ftot = Vector3d();
}

void FixCuttingTool::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixCuttingTool::initial_integrate(Solid &solid, int ip) {
  // cout << "In FixCuttingTool::initial_integrate()\n";
  
  // Not supported
  if (domain->dimension == 3)
    error->one(FLERR, "fix_cuttingtool not supported in 3D\n");

  if (domain->dimension != 2 || !solid.mass.at(ip) || !(solid.mask.at(ip) & groupbit))
    return;

  // Go through all the particles in the group and set b to the right value:
  Vector3d xt(xtvalue .result(mpm, true),
              ytvalue .result(mpm, true),
              ztvalue .result(mpm, true));
  Vector3d vt(vtxvalue.result(mpm, true),
              vtyvalue.result(mpm, true),
              vtzvalue.result(mpm, true));
  Vector3d xA(xAvalue .result(mpm, true),
              yAvalue .result(mpm, true),
              xt.z());
  Vector3d xB(xBvalue .result(mpm, true),
              yBvalue .result(mpm, true),
              xt.z());

  // The equation of line 1 is: (yA - yt) * x - (xA - xt) * y + yt * xA - yA *
  // xt = 0 The equation of line 2 is: (yB - yt) * x - (xB - xt) * y + yt * xB -
  // yB * xt = 0

  double line1[4], line2[4];
  line1[0] = xA[1] - xt[1];
  line1[1] = -xA[0] + xt[0];
  line1[2] = xt[1]*xA[0] - xt[0]*xA[1];
  line1[3] = 1/sqrt(line1[0]*line1[0] + line1[1]*line1[1]);

  line2[0] = xB[1] - xt[1];
  line2[1] = -xB[0] + xt[0];
  line2[2] = xt[1]*xB[0] - xt[0]*xB[1];
  line2[3] = 1/sqrt(line2[0]*line2[0] + line2[1]*line2[1]);

  if (line1[0]*xB[0] + line1[1]*xB[1] + line1[2] < 0)
  {
    line1[0] *= -1;
    line1[1] *= -1;
    line1[2] *= -1;
  }

  if (line2[0]*xA[0] + line2[1]*xA[1] + line2[2] < 0)
  {
    line2[0] *= -1;
    line2[1] *= -1;
    line2[2] *= -1;
  }
  
  Vector3d n1, n2;

  n1[0] = line1[0];
  n1[1] = line1[1];
  n1[2] = 0;
  n1 *= -line1[3];

  n2[0] = line2[0];
  n2[1] = line2[1];
  n2[2] = 0;
  n2 *= -line2[3];

  // cout << "line 1: " << line1[0] << "x + " << line1[1] << "y + " << line1[2] << endl;
  // cout << "line 2: " << line2[0] << "x + " << line2[1] << "y + " << line2[2] << endl;

  double c1p = line1[0]*solid.x.at(ip)[0] + line1[1]*solid.x.at(ip)[1] + line1[1];
  double c2p = line2[0]*solid.x.at(ip)[0] + line2[1]*solid.x.at(ip)[1] + line2[1];

  if (c1p < 0 || c2p < 0)
    return;

  // The particle is inside the tool
  double p1 = fabs(c1p*line1[3]);
  double p2 = fabs(c2p*line2[3]);

  double p;
  Vector3d n;

  if (p1 < p2)
  {
    p = p1;
    n = n1;
  }
  else
  {
    p = p2;
    n = n2;
  }

  const Vector3d &f = K*solid.mat->G*p*(1.0 - solid.damage.at(ip))*n;
  solid.mbp.at(ip) += f;
  ftot += f;
}

void FixCuttingTool::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&K), sizeof(double));
  xtvalue.write_to_restart(of);
  ytvalue.write_to_restart(of);
  ztvalue.write_to_restart(of);

  vtxvalue.write_to_restart(of);
  vtyvalue.write_to_restart(of);
  vtzvalue.write_to_restart(of);

  xAvalue.write_to_restart(of);
  yAvalue.write_to_restart(of);

  xBvalue.write_to_restart(of);
  yBvalue.write_to_restart(of);
}

void FixCuttingTool::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&K), sizeof(double));
  xtvalue.read_from_restart(ifr);
  ytvalue.read_from_restart(ifr);
  ztvalue.read_from_restart(ifr);

  vtxvalue.read_from_restart(ifr);
  vtyvalue.read_from_restart(ifr);
  vtzvalue.read_from_restart(ifr);

  xAvalue.read_from_restart(ifr);
  yAvalue.read_from_restart(ifr);

  xBvalue.read_from_restart(ifr);
  yBvalue.read_from_restart(ifr);
}
