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

#include "output.h"
#include "dump.h"
#include "error.h"
#include "input.h"
#include "log.h"
#include "modify.h"
#include "plot.h"
#include "style_dump.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include "write_restart.h"

#define MIN(A,B) ((A) < (B) ? (A) : (B))

using namespace std;

Output::Output(MPM *mpm) : Pointers(mpm)
{
  ndumps = 0;
  nplots = 0;

  next = 0;
  next_dump_any = 0;
  next_plot_any = 0;

  save_plot = false;

  // create default Log class

  vector<string> log_args;
  log_args.push_back("default");
  log = new Log(mpm, log_args);

  every_log = 1;
  next_log = 1;

  every_restart = next_restart = restart_flag = 0;
  restart = NULL;
}


Output::~Output()
{
  for (int i=0; i<dumps.size();i++) delete dumps[i];
  for (int i=0; i<plots.size();i++) delete plots[i];
  delete log;
  delete restart;
}

void Output::setup(){

  bigint ntimestep = update->ntimestep;

  if (ndumps != 0) {

    if (next_dump.size() != ndumps) next_dump.reserve(ndumps);

    for (int idump = 0; idump < ndumps; idump++) {
      if (every_dump[idump]){
	next_dump[idump] =
          (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      } else {
	cout << "Error every_dump = 0 does not make sense" << endl;
      }

      if (idump != 0) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  if (nplots != 0) {

    if (next_plot.size() != nplots) next_plot.reserve(nplots);
    
    for (int iplot=0; iplot<nplots; iplot++){
      if (every_plot[iplot]){
	next_plot[iplot] =
          (ntimestep/every_plot[iplot])*every_plot[iplot] + every_plot[iplot];
      } else {
	cout << "Error every_plot = 0 does not make sense" << endl;
      }

      if (iplot) next_plot_any = MIN(next_plot_any,next_plot[iplot]);
      else next_plot_any = next_plot[0];
    }
  }

  if (next_plot_any != 0)
    next = MIN(next_dump_any, next_plot_any);
  else
    next = next_dump_any;

  if (restart_flag && next_restart)
    next = MIN(next, next_restart);

  log->init();
  log->header();

  if (every_log) {
    next_log = (ntimestep/every_log)*every_log + every_log;
    if (update->laststep != 0) next_log = MIN(next_log,update->laststep);
  } else next_log = update->laststep;

  if (next!=0) next = MIN(next,next_log);
  else next = next_log;

  if (next==0) {
    error->all(FLERR,"Error: next=0!\n");
  }

  // cout << "next = " << next << endl;
}

void Output::write(bigint ntimestep){

  bigint nsteps = update->nsteps;

  // cout << "In Output::write\n";
  // cout << "restart_flag = " << restart_flag << endl;
  // cout << "next_restart = " << next_restart << endl;
  // cout << "ntimestep = " << ntimestep << endl;

  // If there is at least one dump that requested output at the current step:
  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndumps; idump++) {
      // Which dump requested output:
      if (next_dump[idump] == ntimestep) {
	dumps[idump]->write();
      }

      if (every_dump[idump]) next_dump[idump] += every_dump[idump];

      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  if (restart_flag && next_restart == ntimestep) {
    restart->write();
    next_restart += every_restart;
  }

  if (next_log == ntimestep) {
    modify->run_computes();
    log->write();

    next_log += every_log;
  } else if (ntimestep == 0) {
    modify->run_computes();
    log->write();
  }

  if (restart_flag) {
    next = next_restart;
    if (next_dump_any != 0)
      next = MIN(next, next_dump_any);
    if (next_log != 0)
      next = MIN(next, next_log);
  } else {
    if (next_dump_any != 0)
      next = MIN(next_dump_any, next_log);
    else if (next_log != 0)
      next = next_log;
  }
  
  if (next == 0) {
    error->all(FLERR, "Error: next=0!\n");
  }

  if (next_plot_any == ntimestep) {
    for (int iplot = 0; iplot < nplots; iplot++) {
      // Which plot requested output:
      if (next_plot[iplot] == ntimestep) {
	plots[iplot]->save();
      }

      if (every_plot[iplot]) next_plot[iplot] += every_plot[iplot];

      if (iplot) next_plot_any = MIN(next_plot_any,next_plot[iplot]);
      else next_plot_any = next_plot[0];
    }
  }

  
  if (next_plot_any!=0) next = MIN(next,next_plot_any);
}

void Output::set_log(vector<string> args){
  if (args.size()!=1) {
    error->all(FLERR, "Illegal log command: too many variables.\n");
  }
  every_log = (int) input->parsev(args[0]);
}

void Output::add_dump(vector<string> args){
  cout << "In add_dump" << endl;
  if (args.size() < 5) {
    cout << "Error: not enough arguments in dump command" << endl;
  }

  if (find_dump(args[0]) >= 0) {
    // error->all(FLERR,  "Error: reuse of dump ID.\n");
    delete_dump(args[0]);
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (args[2].compare(#key) == 0) dumps.push_back(new Class(mpm,args));
#include "style_dump.h"
#undef DUMP_CLASS

  else {
    error->all(FLERR, "Unknown dump style " + args[2] + ".\n");
  }

  every_dump.push_back((int) input->parsev(args[3]));
  last_dump.push_back(-1);
  next_dump.push_back(0);
  ndumps++;
}

void Output::modify_dump(vector<string> args){
  cout << "In modify_dump" << endl;

  int idump = find_dump(args[1]);
  if (idump < 0) {
    error->all(FLERR, "Error: dump ID unknown.\n");
  }

  error->all(FLERR, "Unfinished function.\n");
}

void Output::delete_dump(string name){
  cout << "In delete_dump" << endl;

  int idump = find_dump(name);
  if (idump < 0) {
    error->all(FLERR, "Error: dump ID unknown.\n");
  }
  if (dumps[idump]) delete dumps[idump];
  dumps.erase(dumps.begin() + idump);
  every_dump.erase(every_dump.begin() + idump);
  last_dump.erase(last_dump.begin() + idump);
  next_dump.erase(next_dump.begin() + idump);
  ndumps--;
}

int Output::find_dump(string name){
  for (int idump = 0; idump < dumps.size(); idump++){
    if (name.compare(dumps[idump]->id) == 0) return idump;
  }
  return -1;
}

void Output::add_plot(vector<string> args){
  cout << "In add_plot" << endl;
  if (args.size() < 4) {
    cout << "Error: not enough arguments in plot command" << endl;
  }

  if (find_plot(args[0]) >= 0) {
    error->all(FLERR, "Error: reuse of plot ID.\n");
  }

  // create the Plot
  plots.push_back(new Plot(mpm,args));

  every_plot.push_back((int) input->parsev(args[1]));
  last_dump.push_back(-1);
  next_dump.push_back(0);
  nplots++;
}

void Output::modify_plot(vector<string> args){
  cout << "In modify_plot" << endl;

  int iplot = find_plot(args[1]);
  if (iplot == 0) {
    error->all(FLERR, "Error: plot ID unknown.\n");
  }

  error->all(FLERR, "Unfinished function.\n");
}

void Output::delete_plot(string name){
  cout << "In delete_plot" << endl;

  int iplot = find_plot(name);
  if (iplot == 0) {
    error->all(FLERR, "Error: plot ID unknown.\n");
  }

  plots.erase(plots.begin() + iplot);
}

int Output::find_plot(string name){
  for (int iplot = 0; iplot < plots.size(); iplot++){
    if (name.compare(plots[iplot]->id) == 0) return iplot;
  }
  return -1;
}

//  restart(0) => do not write any restart
//  restart(N, file-*.restart) => write a restart every N steps as file-timestep.restart
void Output::create_restart(vector<string> args){
  if (args.size()< 1) {
    error->all(FLERR, "Illegal restart command: too few arguments.\n");
  }
  if (args.size()> 2) {
    error->all(FLERR, "Illegal restart command: too many arguments.\n");
  }

  every_restart = (int) input->parsev(args[0]);
  // cout << "every_restart = " << every_restart << endl;

  if (every_restart == 0) {
    delete restart;
    every_restart = 0;
    next_restart = 0;
    restart_flag = 0;
  } else {
    delete restart;
    restart = new WriteRestart(mpm);
    restart->command(vector<string>(1, args[1]));
    next_restart = (update->ntimestep/every_restart)*every_restart + every_restart;
    // cout << "next_restart = " << next_restart << endl;

    restart_flag = 1;

    if (next) next = MIN(next, next_restart);
  }
}
