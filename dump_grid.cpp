#include <iostream>
#include "output.h"
#include "dump_grid.h"
#include "update.h"
#include "domain.h"
#include "solid.h"
#include "mpmtype.h"

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

    bigint total_nn = 0;
    for (int isolid=0; isolid < domain->solids.size(); isolid++) total_nn += domain->solids[isolid]->grid->nnodes;

    dumpstream << total_nn << endl;
    dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
    dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
    dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
    dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
    dumpstream << "ITEM: ATOMS id type x y z vx vy vz mass\n";

    bigint ID = 0;
    for (int isolid=0; isolid < domain->solids.size(); isolid++) {
      Grid *g = domain->solids[isolid]->grid;
      for (bigint i=0; i<g->nnodes;i++) {
	ID++;
	dumpstream << ID << " " << isolid+1 << " " << g->x[i][0] << " " << g->x[i][1] << " " << g->x[i][2] << " " << g->v[i][0] << " " << g->v[i][1] << " " << g->v[i][2] << " " << g->mass[i] << endl;
      }
    }
    dumpstream.close();
  } else {
    cout << "Error: cannot write in file: " << fdump << endl;
    exit(1);
  }
}
