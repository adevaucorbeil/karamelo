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

#include <dump_particle.h>
#include <domain.h>
#include <error.h>
#include <method.h>
#include <mpm_math.h>
#include <output.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace MPM_Math;

DumpParticle::DumpParticle(MPM *mpm, vector<string> args) : Dump(mpm, args) {
  // cout << "In DumpParticle::DumpParticle()" << endl;
  for (int i = 5; i < args.size(); i++) {
    if (find(known_var.begin(), known_var.end(), args[i]) != known_var.end()) {
      output_var.push_back(args[i]);
    } else {
      string error_str = "Error: output variable \033[1;31m" + args[i] +
                         "\033[0m is unknown!\n";
      error_str += "Availabe output variables: ";
      for (auto v : known_var) {
        error_str += v + ", ";
      }
      error->all(FLERR, error_str);
    }
  }
}

DumpParticle::~DumpParticle()
{
}


void DumpParticle::write()
{
  // Open dump file:
  size_t pos_asterisk = filename.find('*');
  string fdump;

  if (pos_asterisk != string::npos)
    {
      // Replace the asterisk by proc-N.ntimestep:
      fdump = filename.substr(0, pos_asterisk);
      if (universe->nprocs > 1) {
	fdump += "proc-" + to_string(universe->me) + ".";
      }
      fdump += to_string(update->ntimestep);
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
    for (int isolid=0; isolid < domain->solids.size(); isolid++) total_np += domain->solids[isolid]->np_local;

    dumpstream << total_np << endl;
    dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
    dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
    dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
    dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
    dumpstream << "ITEM: ATOMS id type tag ";
    for (auto v: output_var) {
      dumpstream << v << " ";
    }
    dumpstream << endl;

    Matrix3d sigma_;

    for (int isolid=0; isolid < domain->solids.size(); isolid++) {
      Solid *s = domain->solids[isolid];
      for (bigint i=0; i<s->np_local;i++) {
	if (update->method_type == "tlmpm" ||
	    update->method_type == "tlcpdi")
	  sigma_ = s->R[i] * s->sigma[i] * s->R[i].transpose();
	else
	  sigma_ = s->sigma[i];
	dumpstream << s->ptag[i] << " ";
	dumpstream << isolid+1 << " ";
	dumpstream << s->ptag[i] << " ";
	for (auto v: output_var) {
	  if (v.compare("x")==0) dumpstream << s->x[i][0] << " ";
	  else if (v.compare("y")==0) dumpstream << s->x[i][1] << " ";
	  else if (v.compare("z")==0) dumpstream << s->x[i][2] << " ";
	  else if (v.compare("x0")==0) dumpstream << s->x0[i][0] << " ";
	  else if (v.compare("y0")==0) dumpstream << s->x0[i][1] << " ";
	  else if (v.compare("z0")==0) dumpstream << s->x0[i][2] << " ";
	  else if (v.compare("vx")==0) dumpstream << s->v[i][0] << " ";
	  else if (v.compare("vy")==0) dumpstream << s->v[i][1] << " ";
	  else if (v.compare("vz")==0) dumpstream << s->v[i][2] << " ";
	  else if (v.compare("s11")==0) dumpstream << sigma_(0,0) << " ";
	  else if (v.compare("s22")==0) dumpstream << sigma_(1,1) << " ";
	  else if (v.compare("s33")==0) dumpstream << sigma_(2,2) << " ";
	  else if (v.compare("s12")==0) dumpstream << sigma_(0,1) << " ";
	  else if (v.compare("s13")==0) dumpstream << sigma_(0,2) << " ";
	  else if (v.compare("s23")==0) dumpstream << sigma_(1,2) << " ";
	  else if (v.compare("seq")==0) dumpstream << sqrt(3. / 2.) * Deviator(sigma_).norm() << " ";
	  else if (v.compare("e11")==0) dumpstream << s->strain_el[i](0,0) << " ";
	  else if (v.compare("e22")==0) dumpstream << s->strain_el[i](1,1) << " ";
	  else if (v.compare("e33")==0) dumpstream << s->strain_el[i](2,2) << " ";
	  else if (v.compare("e12")==0) dumpstream << s->strain_el[i](0,1) << " ";
	  else if (v.compare("e13")==0) dumpstream << s->strain_el[i](0,2) << " ";
	  else if (v.compare("e23")==0) dumpstream << s->strain_el[i](1,2) << " ";
	  else if (v.compare("damage")==0) dumpstream << s->damage[i] << " ";
	  else if (v.compare("damage_init")==0) dumpstream << s->damage_init[i] << " ";
	  else if (v.compare("volume")==0) dumpstream << s->vol[i] << " ";
	  else if (v.compare("mass")==0) dumpstream << s->mass[i] << " ";
	  else if (v.compare("bx")==0) dumpstream << s->mbp[i][0] << " ";
	  else if (v.compare("by")==0) dumpstream << s->mbp[i][1] << " ";
	  else if (v.compare("bz")==0) dumpstream << s->mbp[i][2] << " ";
	  else if (v.compare("ep")==0) dumpstream << s->eff_plastic_strain[i] << " ";
	  else if (v.compare("epdot")==0) dumpstream << s->eff_plastic_strain_rate[i] << " ";
	  else if (v.compare("ienergy")==0) dumpstream << s->ienergy[i] << " ";
	  else if (v.compare("T")==0) {
	    if (update->method->temp) {
	      dumpstream << s->T[i] << " ";
	    } else {
	      dumpstream << "0 ";
	    }
	  }
	  else if (v.compare("gamma")==0) {
	    if (update->method->temp) {
	      dumpstream << s->gamma[i] << " ";
	    } else {
	      dumpstream << "0 ";	      
	    }
	  }
	}
	dumpstream << endl;
      }
    }
    dumpstream.close();
  } else {
    error->all(FLERR, "Error: cannot write in file: " + fdump + ".\n");
  }
}
