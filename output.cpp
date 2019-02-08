#include "output.h"
#include "dump.h"
#include "input.h"
#include "update.h"
#include "log.h"
#include "style_dump.h"

#define MIN(A,B) ((A) < (B) ? (A) : (B))

using namespace std;

Output::Output(MPM *mpm) : Pointers(mpm)
{
  ndumps = 0;

  // create default Log class

  vector<string> log_args;
  log_args.push_back("default");
  log = new Log(mpm, log_args);

  every_log = 1;
}


Output::~Output()
{
}

void Output::setup(){

  bigint ntimestep = update->ntimestep;

  if (ndumps != 0) {

    if (next_dump.size() != ndumps) next_dump.reserve(ndumps);
    
    for (int idump=0; idump<ndumps; idump++){
      if (every_dump[idump]){
	next_dump[idump] =
          (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      } else {
	cout << "Error every_dump = 0 does not make sense" << endl;
      }

      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  log->init();
  log->header();

  if (every_log) {
    next_log = (ntimestep/every_log)*every_log + every_log;
    next_log = MIN(next_log,update->laststep);
  } else
  next = MIN(next_dump_any,next_log);
}

void Output::write(bigint ntimestep){

  // If there is at least one dump that requested output at the current step:
  if (next_dump_any == ntimestep) {
    for (int idump = 0; idump < ndumps; idump++) {
      // Which dump requested output:
      if (next_dump[idump] == ntimestep) {
	cout << "Should dump .... function to write" << endl;
      }

      if (every_dump[idump]) next_dump[idump] += every_dump[idump];

      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  if (next_log == ntimestep) {
    log->write();

    next_log += every_log;
  }
  
  next = MIN(next_dump_any,next_log);
}

void Output::set_log(vector<string> args){
  if (args.size()!=1) {
    cout << "Illegal log command: too many variables" << endl;
    exit(1);
  }
  every_log = (int) input->parse(args[0]);
}

void Output::add_dump(vector<string> args){
  cout << "In add_dump" << endl;
  if (args.size() < 5) {
    cout << "Error: not enough arguments in dump command" << endl;
  }

  if (find_dump(args[0]) >= 0) {
    cout << "Error: reuse of dump ID" << endl;
    exit(1);
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (args[2].compare(#key) == 0) dumps.push_back(new Class(mpm,args));
#include "style_dump.h"
#undef DUMP_CLASS

  else {
    cout << "Unknown dump style " << args[1] << endl;
    exit(1);
  }

  every_dump.push_back((int) input->parse(args[3]));
  last_dump.push_back(-1);
  ndumps++;
}

void Output::modify_dump(vector<string> args){
  cout << "In modify_dump" << endl;

  int idump = find_dump(args[1]);
  if (idump == 0) {
    cout << "Error: dump ID unknown" << endl;
    exit(1);
  }

  cout << "Unfinished function" << endl;
  exit(1);
}

void Output::delete_dump(string name){
  cout << "In delete_dump" << endl;

  int idump = find_dump(name);
  if (idump == 0) {
    cout << "Error: dump ID unknown" << endl;
    exit(1);
  }

  dumps.erase(dumps.begin() + idump);
}

int Output::find_dump(string name){
  for (int idump = 0; idump < dumps.size(); idump++){
    if (name.compare(dumps[idump]->id) == 0) return idump;
  }
  return -1;
}
