#include <iostream>
#include "output.h"
#include "dump_grid.h"
#include "update.h"
#include "domain.h"
#include "solid.h"
#include "mpmtype.h"
#include <algorithm>

using namespace std;


DumpGrid::DumpGrid(MPM *mpm, vector<string> args) : Dump(mpm, args)
{
  cout << "In DumpGrid::DumpGrid()" << endl;
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
      exit(1);
    }
  }
}

DumpGrid::~DumpGrid()
{
}


void DumpGrid::write()
{
  // Open dump file:
  size_t pos_asterisk = filename.find('*');
  string fdump;

  if (pos_asterisk >= 0)
    {
      // Replace the asterisk by ntimestep:
      fdump = filename.substr(0, pos_asterisk)
	+ to_string(update->ntimestep);
      if (filename.size()-pos_asterisk-1 > 0)
	fdump += filename.substr(pos_asterisk+1, filename.size()-pos_asterisk-1);
    }
  else fdump = filename;

  // cout << "Filemame for dump: " << fdump << endl;

  // Open the file fdump:
  ofstream dumpstream;
  dumpstream.open(fdump, ios_base::out);

  if (dumpstream.is_open()) {
    dumpstream << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n";

    // Check how many different grids we have:
    vector<class Grid *> grids; // We will store the different pointers to grids here.

    for (int isolid=0; isolid < domain->solids.size(); isolid++) {
      if (grids.size()==0) {
	grids.push_back(domain->solids[isolid]->grid);
      } else {
	// If the grid pointer is not present into grids, add it, otherwise continue:
	if ( find(grids.begin(), grids.end(), domain->solids[isolid]->grid) == grids.end() )
	  grids.push_back(domain->solids[isolid]->grid);
      }
    }
    
    // Now loop over the grids to find how many elements there are in total:
    bigint total_nn = 0;
    for (auto g: grids) {
      total_nn += g->nnodes;
    }

    dumpstream << total_nn << endl;
    dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
    dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
    dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
    dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
    dumpstream << "ITEM: ATOMS id type ";
    for (auto v: output_var) {
      dumpstream << v << " ";
    }
    dumpstream << endl;

    bigint ID = 0;
    int igrid = 0;
    for (auto g: grids) {
      for (bigint i=0; i<g->nnodes;i++) {
	ID++;
	dumpstream << i << " "
		   << igrid+1 << " ";
	for (auto v: output_var) {
	  if (v.compare("x")==0) dumpstream << g->x[i][0] << " ";
	  else if (v.compare("y")==0) dumpstream << g->x[i][1] << " ";
	  else if (v.compare("z")==0) dumpstream << g->x[i][2] << " ";
	  else if (v.compare("vx")==0) dumpstream << g->v[i][0] << " ";
	  else if (v.compare("vy")==0) dumpstream << g->v[i][1] << " ";
	  else if (v.compare("vz")==0) dumpstream << g->v[i][2] << " ";
	  else if (v.compare("bx")==0) dumpstream << g->mb[i][0] << " ";
	  else if (v.compare("by")==0) dumpstream << g->mb[i][1] << " ";
	  else if (v.compare("bz")==0) dumpstream << g->mb[i][2] << " ";
	  else if (v.compare("mass")==0) dumpstream << g->mass[i] << " ";
	  else if (v.compare("ntypex")==0) dumpstream << g->ntype[i][0] << " ";
	  else if (v.compare("ntypey")==0) dumpstream << g->ntype[i][1] << " ";
	  else if (v.compare("ntypez")==0) dumpstream << g->ntype[i][2] << " ";
	}
	dumpstream << endl;
      }
    }
    dumpstream.close();
  } else {
    cout << "Error: cannot write in file: " << fdump << endl;
    exit(1);
  }
}
