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

#include <dump_particle_gz.h>
#include <domain.h>
#include <error.h>
#include <method.h>
#include <mpm_math.h>
#include <output.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <gzstream.h>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace MPM_Math;

DumpParticleGz::DumpParticleGz(MPM *mpm, vector<string> args)
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

DumpParticleGz::~DumpParticleGz() {}

void DumpParticleGz::write() {
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

  bigint total_np = 0;
  for (int isolid = 0; isolid < domain->solids.size(); isolid++)
    total_np += domain->solids[isolid]->np_local;

  dumpstream << total_np << endl;
  dumpstream << "ITEM: BOX BOUNDS sm sm sm\n";
  dumpstream << domain->boxlo[0] << " " << domain->boxhi[0] << endl;
  dumpstream << domain->boxlo[1] << " " << domain->boxhi[1] << endl;
  dumpstream << domain->boxlo[2] << " " << domain->boxhi[2] << endl;
  dumpstream << "ITEM: ATOMS id type tag ";
  for (auto v : output_var) {
    dumpstream << v << " ";
  }
  dumpstream << endl;


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
      dumpstream << ptag[isolid][i] << " ";
      dumpstream << isolid + 1 << " ";
      dumpstream << ptag[isolid][i] << " ";
      for (const string &v: output_var)
      {
        if (v == "x")
          dumpstream << x[isolid][i][0] << " ";
        else if (v == "y")
          dumpstream << x[isolid][i][1] << " ";
        else if (v == "z")
          dumpstream << x[isolid][i][2] << " ";
        else if (v == "x0")
          dumpstream << x0[isolid][i][0] << " ";
        else if (v == "y0")
          dumpstream << x0[isolid][i][1] << " ";
        else if (v == "z0")
          dumpstream << x0[isolid][i][2] << " ";
        else if (v == "vx")
          dumpstream << this->v[isolid][i][0] << " ";
        else if (v == "vy")
          dumpstream << this->v[isolid][i][1] << " ";
        else if (v == "vz")
          dumpstream << this->v[isolid][i][2] << " ";
        else if (v == "s11")
          dumpstream << sigma_(0, 0) << " ";
        else if (v == "s22")
          dumpstream << sigma_(1, 1) << " ";
        else if (v == "s33")
          dumpstream << sigma_(2, 2) << " ";
        else if (v == "s12")
          dumpstream << sigma_(0, 1) << " ";
        else if (v == "s13")
          dumpstream << sigma_(0, 2) << " ";
        else if (v == "s23")
          dumpstream << sigma_(1, 2) << " ";
        else if (v == "seq")
          dumpstream << sqrt(3. / 2.) * Deviator(sigma_).norm() << " ";
        else if (v == "e11")
          dumpstream << strain_el[isolid][i](0, 0) << " ";
        else if (v == "e22")
          dumpstream << strain_el[isolid][i](1, 1) << " ";
        else if (v == "e33")
          dumpstream << strain_el[isolid][i](2, 2) << " ";
        else if (v == "e12")
          dumpstream << strain_el[isolid][i](0, 1) << " ";
        else if (v == "e13")
          dumpstream << strain_el[isolid][i](0, 2) << " ";
        else if (v == "e23")
          dumpstream << strain_el[isolid][i](1, 2) << " ";
        else if (v == "damage")
          dumpstream << damage[isolid][i] << " ";
        else if (v == "damage_init")
          dumpstream << damage_init[isolid][i] << " ";
        else if (v == "volume")
          dumpstream << vol[isolid][i] << " ";
        else if (v == "mass")
          dumpstream << mass[isolid][i] << " ";
        else if (v == "bx")
          dumpstream << mbp[isolid][i][0] << " ";
        else if (v == "by")
          dumpstream << mbp[isolid][i][1] << " ";
        else if (v == "bz")
          dumpstream << mbp[isolid][i][2] << " ";
        else if (v == "ep")
          dumpstream << eff_plastic_strain[isolid][i] << " ";
        else if (v == "epdot")
          dumpstream << eff_plastic_strain_rate[isolid][i] << " ";
        else if (v == "ienergy")
          dumpstream << ienergy[isolid][i] << " ";
        else if (v == "T") {
          if (update->method->temp) {
            dumpstream << T[isolid][i] << " ";
          } else {
            dumpstream << "0 ";
          }
        } else if (v == "gamma") {
          if (update->method->temp) {
            dumpstream << gamma[isolid][i] << " ";
          } else {
            dumpstream << "0 ";
          }
        }
      }
      dumpstream << endl;
    }
  }
  dumpstream.close();
}
