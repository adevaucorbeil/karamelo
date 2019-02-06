#include "output.h"
#include "dump.h"
#include "input.h"
#include "style_dump.h"


using namespace std;

Output::Output(MPM *mpm) : Pointers(mpm)
{
  next_dump = NULL;
  last_dump = NULL;
  var_dump = NULL;
  ivar_dump = NULL;
}


Output::~Output()
{
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
