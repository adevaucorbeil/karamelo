#include <iostream>
#include "group.h"
#include "translate_particles.h"
#include "domain.h"
#include "input.h"
#include "var.h"

using namespace std;

TranslateParticles::TranslateParticles(MPM *mpm) : Pointers(mpm) {}

Var TranslateParticles::command(vector<string> args) {
  cout << "In TranslateParticles::command()" << endl;

  if (args.size() < 6) {
    cout << "Illegal translate_particles command. Must be of the form: translate(solid-ID, region, region-ID, delx, dely, delz)" << endl;
    exit(1);
  }

  int ns = domain->solids.size();

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) {
      cout << "Error: solid " << args[0] << " unknown.\n";
      exit(1);
    }
  }

  if (args[1].compare("region")==0) translate_region(args, isolid);
  else {
    cout << "Error: use of illegal keyword for translate_particles command: " << args[1] << endl;
    exit(1);
  }

  return Var(0);
}

void TranslateParticles::translate_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0) {
    cout << "Error: region " << args[2] << " unknown.\n";
    exit(1);
  }

  delx = input->parsev(args[3]);
  if (delx.is_constant() && abs(delx.result())<= 1.0e-12) xset = false;
  else xset = true;
  dely = input->parsev(args[4]);
  if (dely.is_constant() && abs(dely.result())<= 1.0e-12) yset = false;
  else yset = true;
  delz = input->parsev(args[5]);
  if (delz.is_constant() && abs(delz.result())<= 1.0e-12) zset = false;
  else zset = true;

  int ns = domain->solids.size();
  int np;
  double delx_value, dely_value, delz_value;
  Eigen::Vector3d *x0, *x;

  for(int is=0; is<ns; is++){
    if ((isolid < 0) || (is == isolid)) {
      np = domain->solids[is]->np;
      x0 = domain->solids[is]->grid->x0;
      x = domain->solids[is]->grid->x;

      for(int ip=0; ip<np; ip++){
	if (domain->regions[iregion]->inside(x[ip][0], x[ip][1], x[ip][2])==1) {
	  (*input->vars)["x0"] = Var("x0", x0[ip][0]);
	  (*input->vars)["x"] = Var("x", x[ip][0]);
	  (*input->vars)["y0"] = Var("y0", x0[ip][1]);
	  (*input->vars)["y"] = Var("y", x[ip][1]);
	  (*input->vars)["z0"] = Var("z0", x0[ip][2]);
	  (*input->vars)["z"] = Var("z", x[ip][2]);

	  if (xset) {
	    delx_value = delx.result(mpm);
	    x0[ip][0] += delx_value;
	    x[ip][0] += delx_value;
	  }
	  if (yset) {
	    dely_value = dely.result(mpm);
	    x0[ip][1] += dely_value;
	    x[ip][1] += dely_value;
	  }
	  if (zset) {
	    delz_value = delz.result(mpm);
	    x0[ip][2] += delz_value;
	    x[ip][2] += delz_value;
	  }
	}
      }  
    }
  }
}
