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

#include <fix_meld_tool.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace FixConst;


FixMeldTool::FixMeldTool(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];
    return;
  }
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_meldtool needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixMeldTool with ID: " << args[0] << endl;
  }
  id = args[0];

  int k = 2;
  
  K = input->parsev(args[++k]).result(mpm);
  k++;
  if (args[k].compare("x") == 0) {
    dim = X;
    axis0 = Y;
    axis1 = Z;
  } else if (args[k].compare("y") == 0) {
    dim = Y;
    axis0 = X;
    axis1 = Z;
  } else if (args[k].compare("z") == 0) {
    dim = Z;
    axis0 = X;
    axis1 = Y;
  } else {
    error->all(FLERR, "Unknown dim: " + args[k] + "! Options are : x, y, and z.\n");
  }

  w = input->parsev(args[++k]).result(mpm);
  c1 = input->parsev(args[++k]);
  c2 = input->parsev(args[++k]);
  theta = input->parsev(args[++k]);
  lo = input->parsev(args[++k]).result(mpm);
  hi = input->parsev(args[++k]).result(mpm);
  Rmax = input->parsev(args[++k]).result(mpm);
  RmaxSq = Rmax * Rmax;
}

FixMeldTool::~FixMeldTool() {}

void FixMeldTool::init() {}

void FixMeldTool::setup() {}

void FixMeldTool::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixMeldTool::initial_integrate() {
  // cout << "In FixMeldTool::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Vector3d ftot, ftot_reduced;


  double theta_ = theta.result(mpm);

  Vector3d xprime, c;
  Matrix3d R, Rt;

  if (dim == X) {
    R(0,0) = 1;
    R(0,1) = 0;
    R(0,2) = 0;
    R(1,0) = 0;
    R(1,1) = cos(theta_);
    R(1,2) = sin(theta_);
    R(2,0) = 0;
    R(2,1) = -sin(theta_);
    R(2,2) = cos(theta_);

    c(0) = lo;
    c(1) = c1.result(mpm);
    c(2) = c2.result(mpm);
  } else if (dim == Y) {
    R(0,0) = cos(theta_);
    R(0,1) = 0;
    R(0,2) = sin(theta_);
    R(1,0) = 0;
    R(1,1) = 1;
    R(1,2) = 0;
    R(2,0) = -sin(theta_);
    R(2,1) = 0;
    R(2,2) = cos(theta_);

    c(0) = c1.result(mpm);
    c(1) = lo;
    c(2) = c2.result(mpm);
  } else if (dim == Z) {
    R(0,0) = 1;
    R(0,1) = 0;
    R(0,2) = 0;
    R(1,0) = 0;
    R(1,1) = cos(theta_);
    R(1,2) = sin(theta_);
    R(2,0) = 0;
    R(2,1) = -sin(theta_);
    R(2,2) = cos(theta_);

    c(0) = c1.result(mpm);
    c(1) = c2.result(mpm);
    c(2) = lo;
  }

  Rt = R.transpose();
  ftot = Vector3d();

  double p0, p1, p2, pext, p, rSq;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      double KG = K * s->mat->G;

      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mass[ip] > 0) {
	  if (s->mask[ip] & groupbit) {

	    xprime = s->x[ip] - c;
	    // if (update->ntimestep > 89835 && (s->ptag[ip]==12 || s->ptag[ip]==21)) {
	    //   cout << "Check Particle " << s->ptag[ip] << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\n";
	    //   cout << "R=\n" << R << endl;
	    // }
	    if (xprime[dim] >= 0 && xprime[dim] <= hi - lo &&
		xprime(axis0) <= Rmax && xprime(axis0) >= -Rmax &&
		xprime(axis1) <= Rmax && xprime(axis1) >= -Rmax) {
	      // if (s->ptag[ip]==12 || s->ptag[ip]==21) {
	      // 	cout << "Particle " << s->ptag[ip] << " in 1\n";
	      // }

	      rSq = xprime(axis0) * xprime(axis0) + xprime(axis1) * xprime(axis1);

	      if (rSq <= RmaxSq) {
		xprime = R * xprime;
		p0 = xprime[axis0];
		p1 = xprime[axis1];
		p2 = xprime[dim];
		pext = Rmax - sqrt(p0*p0 + p1*p1);

		f = Vector3d();

		if (p0 > w) {
		  p = p0 - w;
		  f[axis0] = -p;
		}

		if (p0 < -w) {
		  p = -w - p0;
		  f[axis0] = p;
		}

		if (p1 > w) {
		  p = p1 - w;
		  f[axis1] = -p;
		}

		if (p1 < -w) {
		  p = -p1 - w;
		  f[axis1] = p;
		}

		if (pext > 0) {
		  if (pext < f.norm()) {
		    f = Vector3d();
		    double r = sqrt(rSq);
		    f[axis0] = pext * xprime[axis0]/r;
		    f[axis1] = pext * xprime[axis1]/r;
		  }
		}

		if (p2 > 0) {
		  if (p2 < f.norm()) {
		    f = Vector3d();
		    f[dim] = -p2;
		  }
		}

		f = KG * (1.0 - s->damage[ip]) * (Rt * f);

		// if (s->ptag[ip]==12 || s->ptag[ip]==21) {
		//   Vector3d dx = s->x[ip] - c;
		//   cout << "Particle " << s->ptag[ip] << " f=[" << f[0] << "," << f[1] << "," << f[2] << "]\tw=" << w << " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\tdx=[" << dx(0) << "," << dx(1) << "," << dx(2) << "]\n";
		//   cout << "R=\n" << R << endl;
		// }
		s->mbp[ip] += f;
		ftot += f;
	      }
	    }
          }
        }
      }
    }
  } else {
    s = domain->solids[solid];
    double KG = K * s->mat->G;

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mass[ip] > 0) {
	if (s->mask[ip] & groupbit) {

	  xprime = s->x[ip] - c;
	  // if (update->ntimestep > 89835 && (s->ptag[ip]==12 || s->ptag[ip]==21)) {
	  //   cout << "Check Particle " << s->ptag[ip] << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\n";
	  //   cout << "R=\n" << R << endl;
	  // }
	  if (xprime[dim] >= 0 && xprime[dim] <= hi-lo &&
	      xprime(axis0) <= Rmax && xprime(axis0) >= -Rmax &&
	      xprime(axis1) <= Rmax && xprime(axis1) >= -Rmax) {
	    // if (s->ptag[ip]==12 || s->ptag[ip]==21) {
	    //   cout << "Particle " << s->ptag[ip] << " in 1\n";
	    // }

	    rSq = xprime(axis0) * xprime(axis0) + xprime(axis1) * xprime(axis1);

	    if (rSq <= RmaxSq) {
	      // if (s->ptag[ip]==12 || s->ptag[ip]==21) {
	      // 	cout << "Particle " << s->ptag[ip] << " in 2\n";
	      // }
	      xprime = R * xprime;
	      p0 = xprime[axis0];
	      p1 = xprime[axis1];
	      p2 = xprime[dim];
		  pext = Rmax - sqrt(p0*p0 + p1*p1);

	      f = Vector3d();

	      if (p0 > w) {
		p = p0 - w;
		f[axis0] = -p;
	      }

	      if (p0 < -w) {
		p = -w - p0;
		f[axis0] = p;
	      }

	      if (p1 > w) {
		p = p1 - w;
		f[axis1] = -p;
	      }

	      if (p1 < -w) {
		p = -p1 - w;
		f[axis1] = p;
	      }

	      if (pext > 0) {
		if (pext < f.norm()) {
		  f = Vector3d();
		  double r = sqrt(rSq);
		  f[axis0] = pext * xprime[axis0]/r;
		  f[axis1] = pext * xprime[axis1]/r;
		}
	      }

	      if (p2 > 0) {
		if (p2 < f.norm()) {
		  f = Vector3d();
		  f[dim] = -p2;
		}
	      }

	      f = KG * (1.0 - s->damage[ip]) * (Rt * f);
	      s->mbp[ip] += f;
	      // if (s->ptag[ip]==12 || s->ptag[ip]==21) {
	      // 	Vector3d dx = s->x[ip] - c;
	      // 	cout << "Particle " << s->ptag[ip] << " f=[" << f[0] << "," << f[1] << "," << f[2] << "]\tw=" << w << " p0=" << p0 << " p1=" << p1 << " p2=" << p2 << "\txprime=[" << xprime[0] << "," << xprime[1] << "," << xprime[2] << "]\tdx=[" << dx(0) << "," << dx(1) << "," << dx(2) << "]\n";
	      // 	cout << "R=\n" << R << endl;
	      // }
	      ftot += f;
	      // if (f[dim] != 0) {
	      // 	cout << "particle " << s->ptag[ip] << " force:" << f[0] << ", " << f[1] << ", " << f[2] << endl;
	      // }
	    }
	  }
	}
      }
    }
  }
  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}


void FixMeldTool::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&dim), sizeof(int));
  of->write(reinterpret_cast<const char *>(&axis0), sizeof(int));
  of->write(reinterpret_cast<const char *>(&axis1), sizeof(int));
  of->write(reinterpret_cast<const char *>(&K), sizeof(double));
  of->write(reinterpret_cast<const char *>(&w), sizeof(double));
  of->write(reinterpret_cast<const char *>(&lo), sizeof(double));
  of->write(reinterpret_cast<const char *>(&hi), sizeof(double));
  of->write(reinterpret_cast<const char *>(&Rmax), sizeof(double));
  
  c1.write_to_restart(of);
  c2.write_to_restart(of);
  theta.write_to_restart(of);
}

void FixMeldTool::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&dim), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&axis0), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&axis1), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&K), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&w), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&lo), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&hi), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&Rmax), sizeof(double));
  RmaxSq = Rmax * Rmax;

  c1.read_from_restart(ifr);
  c2.read_from_restart(ifr);
  theta.read_from_restart(ifr);
}
