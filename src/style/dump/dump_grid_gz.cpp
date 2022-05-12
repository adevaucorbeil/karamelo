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

#include <dump_grid_gz.h>
#include <domain.h>
#include <error.h>
#include <output.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <algorithm>
#include <gzstream.h>
#include <iostream>

using namespace std;

DumpGridGz::DumpGridGz(MPM *mpm, vector<string> args) : Dump(mpm, args) {
  cout << "In DumpGridGz::DumpGridGz()" << endl;
  for (int i=5; i<args.size(); i++){
    if (find(known_var.begin(), known_var.end(), args[i]) != known_var.end()) {
      output_var.push_back(args[i]);
    } else {
      cout << "Error: output variable \033[1;31m" << args[i] << "\033[0m is unknown!\n";
      cout << "Availabe output variables: ";
      for (auto v: known_var) {
	cout << v << ", ";
      }
      cout << endl;
      error->all(FLERR, "");
    }
  }

  bool xyz = true;
  bool vxyz = true;
  bool bxyz = true;
  bool ntypexyz = true;

  const Grid &g = *domain->grid;

  ntag = create_mirror(g.ntag);

  for (const string &v: output_var)
  {
    if ((v == "x" || v == "y" || v == "z") && xyz)
    {
      x = create_mirror(g.x);
      xyz = false;
      cout << "x.extent=" << x.extent(0) << endl;
    }
    else if ((v == "vx" || v == "vy" || v == "vz") && vxyz)
    {
      this->v = create_mirror(g.v);
      vxyz = false;
    }
    else if ((v == "bx" || v == "by" || v == "bz") && bxyz)
    {
      mb = create_mirror(g.mb);
      bxyz = false;
    }
    else if (v == "mass")
      mass = create_mirror(g.mass);
    else if (v.compare("mask")==0)
      mask = create_mirror(g.mask);
    else if ((v == "ntypex" || v == "ntypey" || v == "ntypez") && ntypexyz)
    {
      ntype = create_mirror(g.ntype);
      ntypexyz = false;
    }
    else if (v == "T")
      T = create_mirror(g.T);
    else if (v == "rigid")
      rigid = create_mirror(g.rigid);
  }
}

DumpGridGz::~DumpGridGz() {}

void DumpGridGz::write() {
  // Open dump file:
  size_t pos_asterisk = filename.find('*');
  string fdump;

  if (pos_asterisk != string::npos) {
    // Replace the asterisk by proc-N.ntimestep:
    fdump = filename.substr(0, pos_asterisk);
    if (universe->nprocs > 1) {
      fdump += "proc-" + to_string(universe->me) + ".";
    }
    fdump += to_string(update->ntimestep);
    if (filename.size() - pos_asterisk - 1 > 0)
      fdump +=
          filename.substr(pos_asterisk + 1, filename.size() - pos_asterisk - 1);
  } else
    fdump = filename;

  // cout << "Filemame for dump: " << fdump << endl;

  // Open the file fdump:
  ogzstream dumpstream(fdump.c_str());

  dumpstream << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n";

  // Check how many different grids we have:
  vector<class Grid *>
      grids; // We will store the different pointers to grids here.

  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    if (grids.size() == 0) {
      grids.push_back(domain->solids[isolid]->grid);
    } else {
      // If the grid pointer is not present into grids, add it, otherwise
      // continue:
      if (find(grids.begin(), grids.end(), domain->solids[isolid]->grid) ==
          grids.end())
        grids.push_back(domain->solids[isolid]->grid);
    }
  }

  // Now loop over the grids to find how many elements there are in total:
  bigint total_nn = 0;
  for (auto g : grids) {
    total_nn += g->nnodes_local;
  }

  dumpstream << total_nn << endl;
  dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
  dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
  dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
  dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
  dumpstream << "ITEM: ATOMS id type ";
  for (auto v : output_var) {
    dumpstream << v << " ";
  }
  dumpstream << endl;

  int igrid = 0;
  for (auto g: grids) {
    bool xyz = true;
    bool vxyz = true;
    bool bxyz = true;
    bool ntypexyz = true;

    deep_copy(ntag, g->ntag);

    for (const string &v: output_var)
    {
      if ((v == "x" || v == "y" || v == "z") && xyz)
      {
        deep_copy(x, g->x);
        xyz = false;
      }
      else if ((v == "vx" || v == "vy" || v == "vz") && vxyz)
      {
        deep_copy(this->v, g->v);
        vxyz = false;
      }
      else if ((v == "bx" || v == "by" || v == "bz") && bxyz)
      {
        deep_copy(mb, g->mb);
        bxyz = false;
      }
      else if (v == "mass")
        deep_copy(mass, g->mass);
      else if (v.compare("mask")==0)
        deep_copy(mask, g->mask);
      else if ((v == "ntypex" || v == "ntypey" || v == "ntypez") && ntypexyz)
      {
        deep_copy(ntype, g->ntype);
        ntypexyz = false;
      }
      else if (v == "T")
        deep_copy(T, g->T);
      else if (v == "rigid")
        deep_copy(rigid, g->rigid);
    }

    for (bigint i = 0; i < g->nnodes_local; i++) {
      dumpstream << ntag[i] << " " << igrid + 1 << " ";
      for (auto v : output_var) {
        if (v == "x")
          dumpstream << x[i][0] << " ";
        else if (v == "y")
          dumpstream << x[i][1] << " ";
        else if (v == "z")
          dumpstream << x[i][2] << " ";
        else if (v == "vx")
          dumpstream << this->v[i][0] << " ";
        else if (v == "vy")
          dumpstream << this->v[i][1] << " ";
        else if (v == "vz")
          dumpstream << this->v[i][2] << " ";
        else if (v == "bx")
          dumpstream << mb[i][0] << " ";
        else if (v == "by")
          dumpstream << mb[i][1] << " ";
        else if (v == "bz")
          dumpstream << mb[i][2] << " ";
        else if (v == "mass")
          dumpstream << mass[i] << " ";
	      else if (v.compare("mask")==0)
	        dumpstream << mask[i] << " ";
        else if (v == "ntypex")
          dumpstream << ntype[i][0] << " ";
        else if (v == "ntypey")
          dumpstream << ntype[i][1] << " ";
        else if (v == "ntypez")
          dumpstream << ntype[i][2] << " ";
        else if (v == "T")
          dumpstream << T[i] << " ";
        else if (v == "rigid")
          dumpstream << rigid[i] << " ";
      }
      dumpstream << endl;
    }
  }
  dumpstream.close();
}
