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
  // cout << "In DumpGridGz::DumpGridGz()" << endl;
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
  for (auto g : grids) {
    for (bigint i = 0; i < g->nnodes_local; i++) {
      dumpstream << g->ntag[i] << " " << igrid + 1 << " ";
      for (auto v : output_var) {
        if (v == "x")
          dumpstream << g->x[i][0] << " ";
        else if (v == "y")
          dumpstream << g->x[i][1] << " ";
        else if (v == "z")
          dumpstream << g->x[i][2] << " ";
        else if (v == "vx")
          dumpstream << g->v[i][0] << " ";
        else if (v == "vy")
          dumpstream << g->v[i][1] << " ";
        else if (v == "vz")
          dumpstream << g->v[i][2] << " ";
        else if (v == "bx")
          dumpstream << g->mb[i][0] << " ";
        else if (v == "by")
          dumpstream << g->mb[i][1] << " ";
        else if (v == "bz")
          dumpstream << g->mb[i][2] << " ";
        else if (v == "mass")
          dumpstream << g->mass[i] << " ";
	else if (v.compare("mask")==0)
	  dumpstream << g->mask[i] << " ";
        else if (v == "ntypex")
          dumpstream << g->ntype[i][0] << " ";
        else if (v == "ntypey")
          dumpstream << g->ntype[i][1] << " ";
        else if (v == "ntypez")
          dumpstream << g->ntype[i][2] << " ";
        else if (v == "T")
          dumpstream << g->T[i] << " ";
        else if (v == "rigid")
          dumpstream << g->rigid[i] << " ";
      }
      dumpstream << endl;
    }
  }
  dumpstream.close();
}
