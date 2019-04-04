#include <iostream>
#include "output.h"
#include "dump_particle.h"
#include "update.h"
#include "domain.h"
#include "solid.h"
#include "mpmtype.h"
#include "mpm_math.h"

using namespace std;
using namespace MPM_Math;


DumpParticle::DumpParticle(MPM *mpm, vector<string> args) : Dump(mpm, args)
{
  cout << "In DumpParticle::DumpParticle()" << endl;
}

DumpParticle::~DumpParticle()
{
}


void DumpParticle::write()
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
  ofstream dumpstream(fdump, ios_base::out);

  if (dumpstream.is_open()) {
    dumpstream << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n";

    bigint total_np = 0;
    for (int isolid=0; isolid < domain->solids.size(); isolid++) total_np += domain->solids[isolid]->np;

    dumpstream << total_np << endl;
    dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
    dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
    dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
    dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
    dumpstream << "ITEM: ATOMS id type x y z vx vy vz s11 s22 s33 s12 s13 s23 seq damage damage_init volume mass\n";

    bigint ID = 0;
    for (int isolid=0; isolid < domain->solids.size(); isolid++) {
      Solid *s = domain->solids[isolid];
      for (bigint i=0; i<s->np;i++) {
	ID++;
	dumpstream << ID << " ";
	dumpstream << isolid+1 << " ";
	dumpstream << s->x[i][0] << " ";
	dumpstream << s->x[i][1] << " ";
	dumpstream << s->x[i][2] << " ";
	dumpstream << s->v[i][0] << " ";
	dumpstream << s->v[i][1] << " ";
	dumpstream << s->v[i][2] << " ";
	dumpstream << s->sigma[i](0,0) << " ";
	dumpstream << s->sigma[i](1,1) << " ";
	dumpstream << s->sigma[i](2,2) << " ";
	dumpstream << s->sigma[i](0,1) << " ";
	dumpstream << s->sigma[i](0,2) << " ";
	dumpstream << s->sigma[i](1,2) << " ";
	dumpstream << sqrt(3. / 2.) * Deviator(s->sigma[i]).norm() << " ";
	dumpstream << s->damage[i] << " ";
	dumpstream << s->damage_init[i] << " ";
	dumpstream << s->vol[i] << " ";
	dumpstream << s->mass[i] << endl;
      }
    }
    dumpstream.close();
  } else {
    cout << "Error: cannot write in file: " << fdump << endl;
    exit(1);
  }
}
