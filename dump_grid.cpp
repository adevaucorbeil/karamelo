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
    dumpstream << "ITEM: ATOMS id type x y z vx vy vz fbx fby fbz mass ntypex ntypey ntypez\n";

    bigint ID = 0;
    int igrid = 0;
    for (auto g: grids) {
      for (bigint i=0; i<g->nnodes;i++) {
	ID++;
	dumpstream << ID << " "
		   << igrid+1 << " "
		   << g->x[i][0] << " " << g->x[i][1] << " " << g->x[i][2] << " "
		   << g->v[i][0] << " " << g->v[i][1] << " " << g->v[i][2] << " "
		   << g->b[i][0] << " " << g->b[i][1] << " " << g->b[i][2] << " "
		   << g->mass[i] << " "
		   << g->ntype[i][0] << " " << g->ntype[i][1] << " " << g->ntype[i][2] << endl;
      }
    }
    dumpstream.close();
  } else {
    cout << "Error: cannot write in file: " << fdump << endl;
    exit(1);
  }
}
