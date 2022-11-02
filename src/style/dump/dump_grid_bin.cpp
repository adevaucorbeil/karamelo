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

#include <dump_grid_bin.h>
#include <domain.h>
#include <error.h>
#include <output.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <algorithm>
#include <iostream>

using namespace std;

DumpGridBin::DumpGridBin(MPM *mpm, vector<string> args) : Dump(mpm, args) {
  // cout << "In DumpGridBin::DumpGridBin()" << endl;
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
  bool vxyz_update = true;
  bool bxyz = true;
  bool ntypexyz = true;

  const Grid &g  = *domain->solids.front()->grid;

  ntag = create_mirror(g.ntag);

  for (const string &v: output_var)
  {
    if ((v == "x" || v == "y" || v == "z") && xyz)
    {
      x = create_mirror(g.x);
      xyz = false;
    }
    else if ((v == "vx" || v == "vy" || v == "vz") && vxyz)
    {
      this->v = create_mirror(g.v);
      vxyz = false;
    }
    else if ((v == "vx_update" || v == "vy_update" || v == "vz_update") && vxyz_update)
    {
      this->v_update = create_mirror(g.v_update);
      vxyz_update = false;
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
    else if (v.compare("vol")==0)
      vol = create_mirror(g.vol);
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

DumpGridBin::~DumpGridBin() {
  // Wait for all threads to be completed.
  for (int i=0; i<threads.size(); i++)
    threads[i].first.join();
}

void DumpGridBin::write() {

  int ithread;
  pair<thread, vector<float>> *th = nullptr;

  for (int i=0; i<threads.size(); i++) {
    if (threads[i].second.empty()) {
      th = &threads[i];
      threads[i].first.join();
      ithread = i;
      break;
    }
  }

  if (!th) {
    ithread = threads.size();
    threads.emplace_back();
    th = &threads.back();
  }

  vector<float> &buf = th->second;

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

  // Check how many different grids we have:
  vector<class Grid *>  grids; // We will store the different pointers to grids here.

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
    total_nn += g->nnodes_local * g->nsolids;
  }

  int size_one = output_var.size() + 2;

  // Resize the buffer:
  buf.reserve(total_nn * size_one);

  int igrid = 0;
  for (auto g: grids) {
    bool xyz = true;
    bool vxyz = true;
    bool vxyz_update = true;
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
      else if ((v == "vx_update" || v == "vy_update" || v == "vz_update") && vxyz_update)
      {
        deep_copy(this->v_update, g->v_update);
        vxyz_update = false;
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
      else if (v.compare("vol")==0)
        deep_copy(vol, g->vol);
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

    for (int is = 0; is < g->nsolids; is++) {
      for (bigint i = 0; i < g->nnodes_local; i++) {
	buf.push_back(ntag[i]);
	buf.push_back(igrid + 1 + is);
	for (auto v : output_var) {
	  if (v == "x")
	    buf.push_back(x[i][0]);
	  else if (v == "y")
	    buf.push_back(x[i][1]);
	  else if (v == "z")
	    buf.push_back(x[i][2]);
	  else if (v == "vx")
	    buf.push_back(this->v(is, i)[0]);
	  else if (v == "vy")
	    buf.push_back(this->v(is, i)[1]);
	  else if (v == "vz")
	    buf.push_back(this->v(is, i)[2]);
	  else if (v == "vx_update")
	    buf.push_back(this->v_update(is, i)[0]);
	  else if (v == "vy_update")
	    buf.push_back(this->v_update(is, i)[1]);
	  else if (v == "vz_update")
	    buf.push_back(this->v_update(is, i)[2]);
	  else if (v == "bx")
	    buf.push_back(mb(is, i)[0]);
	  else if (v == "by")
	    buf.push_back(mb(is, i)[1]);
	  else if (v == "bz")
	    buf.push_back(mb(is, i)[2]);
	  else if (v == "mass")
	    buf.push_back(mass(is, i));
	  else if (v.compare("mask")==0)
	    buf.push_back(mask[i]);
	  else if (v.compare("vol")==0)
	    buf.push_back(vol(is, i));
	  else if (v == "ntypex")
	    buf.push_back(ntype[i][0]);
	  else if (v == "ntypey")
	    buf.push_back(ntype[i][1]);
	  else if (v == "ntypez")
	    buf.push_back(ntype[i][2]);
	  else if (v == "T")
	    buf.push_back(T(is, i));
	  else if (v == "rigid")
	    buf.push_back(rigid(is, i));
	}
      }
    }
  }

  th->first = thread(&DumpGridBin::write_to_file, this, ithread, fdump, total_nn, update->ntimestep);
}

void DumpGridBin::write_to_file(bigint i, string fdump, bigint total_nn, bigint timestep) {
  // Open the file fdump:
  ofstream dumpstream;
  dumpstream.open(fdump.c_str(), std::fstream::out | std::fstream::binary | std::fstream::trunc);

  if (!dumpstream) {
    error->one(FLERR, "Cannot open file " + fdump + ".\n");
  }

  // use negative ntimestep as marker for new format
  bigint fmtlen = MAGIC_STRING.length();
  bigint marker = -fmtlen;
  dumpstream.write(reinterpret_cast<const char *>(&marker),sizeof(bigint));
  dumpstream.write(reinterpret_cast<const char *>(MAGIC_STRING.c_str()), MAGIC_STRING.size());
  dumpstream.write(reinterpret_cast<const char *>(&ENDIAN),sizeof(int));
  dumpstream.write(reinterpret_cast<const char *>(&FORMAT_REVISION),sizeof(int));


  dumpstream.write(reinterpret_cast<const char *>(&timestep),sizeof(bigint));


  dumpstream.write(reinterpret_cast<const char *>(&total_nn),sizeof(bigint));
  int triclinic = 0;
  dumpstream.write(reinterpret_cast<const char *>(&triclinic),sizeof(int));

  int one = 1;
  for(int i=0; i<6; i++)
    dumpstream.write(reinterpret_cast<const char *>(&one),sizeof(int)); // Boundary types

  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[0]),sizeof(float));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[0]),sizeof(float));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[1]),sizeof(float));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[1]),sizeof(float));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxlo[2]),sizeof(float));
  dumpstream.write(reinterpret_cast<const char *>(&domain->boxhi[2]),sizeof(float));
  int size_one = output_var.size() + 2;
  dumpstream.write(reinterpret_cast<const char *>(&size_one),sizeof(int));


  // We are not setting any unit style, so we write 0:
  int unit_style = 0;
  dumpstream.write(reinterpret_cast<const char *>(&unit_style),sizeof(int));

  // We are not storing the simulation time, so we write 0:
  char time_flag = 0;
  dumpstream.write(reinterpret_cast<const char *>(&time_flag),sizeof(char));

  // Write column names:
  string columns = "id type ";
  for (auto v: output_var) {
    columns +=  v + " ";
  }
  int Nc = strlen(columns.c_str());
  dumpstream.write(reinterpret_cast<const char *>(&Nc), sizeof(int));
  dumpstream.write(columns.data(), Nc*sizeof(char));


  int nclusterprocs = 1;
  dumpstream.write(reinterpret_cast<const char *>(&nclusterprocs),sizeof(int));

  int nme = (int) (total_nn * size_one); // # of dump lines this proc contributes to dump (nn->local since each cpu creates its own dump.
  dumpstream.write(reinterpret_cast<const char *>(&nme),sizeof(int));

  vector<float> &buf = threads[i].second;
  dumpstream.write(reinterpret_cast<const char *>(&buf[0]),buf.size()*sizeof(float));

  dumpstream.close();

  // Empty buffer:
  buf.clear();
}
