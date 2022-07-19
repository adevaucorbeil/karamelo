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

#include <dump_particle_bin.h>
#include <algorithm>
#include <domain.h>
#include <error.h>
#include <fstream>
#include <iostream>
#include <matrix.h>
#include <method.h>
#include <mpm_math.h>
#include <output.h>
#include <solid.h>
#include <universe.h>
#include <update.h>

using namespace std;
using namespace MPM_Math;

DumpParticleBin::DumpParticleBin(MPM *mpm, vector<string> args)
    : Dump(mpm, args) {
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


  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    const Solid &s = *domain->solids[isolid];

    bool xyz = true;
    bool x0yz = true;
    bool vxyz = true;
    bool sxyz = true;
    bool exyz = true;
    bool bxyz = true;

    ptag.emplace_back(create_mirror(s.ptag));

    for (const string &v: output_var)
      {
	if ((v == "x" || v == "y" || v == "z") && xyz)
	  {
	    x.emplace_back(create_mirror(s.x));
	    xyz = false;
	  }
	else if ((v == "x0" || v == "y0" || v == "z0") && x0yz)
	  {
	    x0.emplace_back(create_mirror(s.x0));
	    x0yz = false;
	  }
	else if ((v == "vx" || v == "vy" || v == "vz") && vxyz)
	  {
	    this->v.emplace_back(create_mirror(s.v));
	    vxyz = false;
	  }
	else if ((v == "s11" || v == "s22" || v == "s33" ||
		  v == "s12" || v == "s13" || v == "s23" || v == "seq") && sxyz)
	  {
	    sigma.emplace_back(create_mirror(s.sigma));

	    if (update->method_type == "tlmpm" ||
		update->method_type == "tlcpdi")
	      R.emplace_back(create_mirror(s.R));

	    sxyz = false;
	  }
	else if ((v == "e11" || v == "e22" || v == "e33" ||
		  v == "e12" || v == "e13" || v == "e23") && exyz)
	  {
	    strain_el.emplace_back(create_mirror(s.strain_el));
	    exyz = false;
	  }
	else if (v == "damage")
	  damage.emplace_back(create_mirror(s.damage));
	else if (v == "damage_init")
	  damage_init.emplace_back(create_mirror(s.damage_init));
	else if (v == "volume")
	  vol.emplace_back(create_mirror(s.vol));
	else if (v == "mass")
	  mass.emplace_back(create_mirror(s.mass));
	else if ((v == "bx" || v == "by" || v == "bz") && bxyz)
	  {
	    mbp.emplace_back(create_mirror(s.mbp));
	    bxyz = false;
	  }
	else if (v == "ep")
	  eff_plastic_strain.emplace_back(create_mirror(s.eff_plastic_strain));
	else if (v == "epdot")
	  eff_plastic_strain_rate.emplace_back(create_mirror(s.eff_plastic_strain_rate));
	else if (v == "ienergy")
	  ienergy.emplace_back(create_mirror(s.ienergy));
	else if (v == "T" && update->method->temp)
	  T.emplace_back(create_mirror(s.T));
	else if (v == "gamma" && update->method->temp)
	  gamma.emplace_back(create_mirror(s.gamma));
      }
  }
}

DumpParticleBin::~DumpParticleBin() {}

void DumpParticleBin::write() {
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
  ofstream dumpstream;
  dumpstream.open(fdump.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

  if (!dumpstream) {
    error->one(FLERR, "Cannot open file " + fdump + ".\n");
  }

  dumpstream.write(reinterpret_cast<const char *>(&update->ntimestep),sizeof(bigint));

  bigint total_np = 0;
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
    total_np += domain->solids[isolid]->np_local;

  dumpstream.write(reinterpret_cast<const char *>(&total_np),sizeof(bigint));
  int triclinic = 0;
  dumpstream.write(reinterpret_cast<const char *>(&triclinic),sizeof(int));

  int one = 1;
  for(int i=0; i<6; i++)
    dumpstream.write(reinterpret_cast<const char *>(&one),sizeof(int));

  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[0]),sizeof(double));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[0]),sizeof(double));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[1]),sizeof(double));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[1]),sizeof(double));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[2]),sizeof(double));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[2]),sizeof(double));
  int size_one = output_var.size() + 2;
  dumpstream.write(reinterpret_cast<const char *>(&size_one),sizeof(int));
  int nclusterprocs = 1;
  dumpstream.write(reinterpret_cast<const char *>(&nclusterprocs),sizeof(int));

  int nme = (int) total_np; // # of dump lines this proc contributes to dump (np->local since each cpu creates its own dump.
  dumpstream.write(reinterpret_cast<const char *>(&nme),sizeof(int));

  double buf[total_np*size_one];

  int m = 0;
  Matrix3d sigma_;

  for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
    Solid *s = domain->solids[isolid];

    bool xyz = true;
    bool x0yz = true;
    bool vxyz = true;
    bool sxyz = true;
    bool exyz = true;
    bool bxyz = true;

    deep_copy(ptag[isolid], s->ptag);

    for (const string &v: output_var)
    {
      if ((v == "x" || v == "y" || v == "z") && xyz)
      {
        deep_copy(x[isolid], s->x);
        xyz = false;
      }
      else if ((v == "x0" || v == "y0" || v == "z0") && x0yz)
      {
        deep_copy(x0[isolid], s->x0);
        x0yz = false;
      }
      else if ((v == "vx" || v == "vy" || v == "vz") && vxyz)
      {
        deep_copy(this->v[isolid], s->v);
        vxyz = false;
      }
      else if ((v == "s11" || v == "s22" || v == "s33" ||
                v == "s12" || v == "s13" || v == "s23" || v == "seq") && sxyz)
      {
        deep_copy(sigma[isolid], s->sigma);

        if (update->method_type == "tlmpm" ||
            update->method_type == "tlcpdi")
          deep_copy(R[isolid], s->R);

        sxyz = false;
      }
      else if ((v == "e11" || v == "e22" || v == "e33" ||
                v == "e12" || v == "e13" || v == "e23") && exyz)
      {
        deep_copy(strain_el[isolid], s->strain_el);
        exyz = false;
      }
      else if (v == "damage")
        deep_copy(damage[isolid], s->damage);
      else if (v == "damage_init")
        deep_copy(damage_init[isolid], s->damage_init);
      else if (v == "volume")
        deep_copy(vol[isolid], s->vol);
      else if (v == "mass")
        deep_copy(mass[isolid], s->mass);
      else if ((v == "bx" || v == "by" || v == "bz") && bxyz)
      {
        deep_copy(mbp[isolid], s->mbp);
        bxyz = false;
      }
      else if (v == "ep")
        deep_copy(eff_plastic_strain[isolid], s->eff_plastic_strain);
      else if (v == "epdot")
        deep_copy(eff_plastic_strain_rate[isolid], s->eff_plastic_strain_rate);
      else if (v == "ienergy")
        deep_copy(ienergy[isolid], s->ienergy);
      else if (v == "T" && update->method->temp)
        deep_copy(T[isolid], s->T);
      else if (v == "gamma" && update->method->temp)
        deep_copy(gamma[isolid], s->gamma);
    }

    for (bigint i = 0; i < s->np_local; i++) {
      if (update->method_type == "tlmpm" ||
	        update->method_type == "tlcpdi")
	      sigma_ = R[isolid][i]*sigma[isolid][i]*R[isolid][i].transpose();
      else
	      sigma_ = sigma[isolid][i];

      buf[m++] = (double) ptag[isolid][i];
      buf[m++] = (double) isolid + 1;
      for (const string &v: output_var)
      {
        if (v == "x")
          buf[m++] = (double) x[isolid][i][0];
        else if (v == "y")
          buf[m++] = (double) x[isolid][i][1];
        else if (v == "z")
          buf[m++] = (double) x[isolid][i][2];
        else if (v == "x0")
          buf[m++] = (double) x0[isolid][i][0];
        else if (v == "y0")
          buf[m++] = (double) x0[isolid][i][1];
        else if (v == "z0")
          buf[m++] = (double) x0[isolid][i][2];
        else if (v == "vx")
          buf[m++] = (double) this->v[isolid][i][0];
        else if (v == "vy")
          buf[m++] = (double) this->v[isolid][i][1];
        else if (v == "vz")
          buf[m++] = (double) this->v[isolid][i][2];
        else if (v == "s11")
          buf[m++] = (double) sigma_(0, 0);
        else if (v == "s22")
          buf[m++] = (double) sigma_(1, 1);
        else if (v == "s33")
          buf[m++] = (double) sigma_(2, 2);
        else if (v == "s12")
          buf[m++] = (double) sigma_(0, 1);
        else if (v == "s13")
          buf[m++] = (double) sigma_(0, 2);
        else if (v == "s23")
          buf[m++] = (double) sigma_(1, 2);
        else if (v == "seq")
          buf[m++] = (double) sqrt(3. / 2.) * Deviator(sigma_).norm();
        else if (v == "e11")
          buf[m++] = (double) strain_el[isolid][i](0, 0);
        else if (v == "e22")
          buf[m++] = (double) strain_el[isolid][i](1, 1);
        else if (v == "e33")
          buf[m++] = (double) strain_el[isolid][i](2, 2);
        else if (v == "e12")
          buf[m++] = (double) strain_el[isolid][i](0, 1);
        else if (v == "e13")
          buf[m++] = (double) strain_el[isolid][i](0, 2);
        else if (v == "e23")
          buf[m++] = (double) strain_el[isolid][i](1, 2);
        else if (v == "damage")
          buf[m++] = (double) damage[isolid][i];
        else if (v == "damage_init")
          buf[m++] = (double) damage_init[isolid][i];
        else if (v == "volume")
          buf[m++] = (double) vol[isolid][i];
        else if (v == "mass")
          buf[m++] = (double) mass[isolid][i];
        else if (v == "bx")
          buf[m++] = (double) mbp[isolid][i][0];
        else if (v == "by")
          buf[m++] = (double) mbp[isolid][i][1];
        else if (v == "bz")
          buf[m++] = (double) mbp[isolid][i][2];
        else if (v == "ep")
          buf[m++] = (double) eff_plastic_strain[isolid][i];
        else if (v == "epdot")
          buf[m++] = (double) eff_plastic_strain_rate[isolid][i];
        else if (v == "ienergy")
          buf[m++] = (double) ienergy[isolid][i];
        else if (v == "T") {
          if (update->method->temp) {
            buf[m++] = (double) T[isolid][i];
          } else {
            buf[m++] = (double) 0.0;
          }
        } else if (v == "gamma") {
          if (update->method->temp) {
            buf[m++] = (double) gamma[isolid][i];
          } else {
            buf[m++] = (double) 0.0;
          }
        }
      }
    }
  }

  if (m != total_np*size_one) {
    error->one(FLERR, "m != total_np*size_one\n");
  }

  dumpstream.write(reinterpret_cast<const char *>(&buf[0]),m*sizeof(double));
  
  dumpstream.close();
}
