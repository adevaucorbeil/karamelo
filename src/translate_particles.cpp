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

#include "translate_particles.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "var.h"
#include <iostream>

using namespace std;

TranslateParticles::TranslateParticles(MPM *mpm) : Pointers(mpm) {}

Var TranslateParticles::command(vector<string> args) {
  // cout << "In TranslateParticles::command()" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  int ns = domain->solids.size();

  int isolid = domain->find_solid(args[0]);

  if (isolid < 0) {
    if (args[0].compare("all")!=0) {
      cout << "Error: solid " << args[0] << " unknown.\n";
      error->all(FLERR, "");
    }
  }

  if (args[1].compare("region")==0) translate_region(args, isolid);
  else {
    cout << "Error: use of illegal keyword for translate_particles command: " << args[1] << endl;
    error->all(FLERR, "");
  }

  return Var(0);
}

void TranslateParticles::translate_region(vector<string> args, int isolid) {
  int iregion = domain->find_region(args[2]);

  if (iregion < 0)
    {
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
  Solid *s;
  Grid *g;

  for(int is = 0; is < ns; is++)
    {
      if ((isolid < 0) || (is == isolid))
	{
	  s = domain->solids[is];
	  g = s->grid;

	  for(int ip=0; ip < s->np_local; ip++)
	    {
	      if (domain->regions[iregion]->inside(s->x[ip][0], s->x[ip][1], s->x[ip][2])==1)
		{
		  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
		  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
		  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

		  (*input->vars)["x"]  = Var("x", s->x[ip][0]);
		  (*input->vars)["y"]  = Var("y", s->x[ip][1]);
		  (*input->vars)["z"]  = Var("z", s->x[ip][2]);

		  if (xset)
		    {
		      delx_value    = delx.result(mpm);
		      s->x0[ip][0] += delx_value;
		      s->x[ip][0]  += delx_value;
		    }
		  if (yset)
		    {
		      dely_value    = dely.result(mpm);
		      s->x0[ip][1] += dely_value;
		      s->x[ip][1]  += dely_value;
		    }
		  if (zset)
		    {
		      delz_value    = delz.result(mpm);
		      s->x0[ip][2] += delz_value;
		      s->x[ip][2]  += delz_value;
		    }
		}
	    }  
	}
    }
}
