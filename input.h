/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_INPUT_H
#define MPM_INPUT_H

#include "pointers.h"
#include <vector>
#include <map>

class Input : protected Pointers {
public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command

  map<string, class Var> *vars; // List of global variables.

  // functions
  Input(class MPM *, int, char **);
  ~Input();
  void file();                 // process all input
  class Var parsev(string);    // parse an input text line
  double parse(string);        // deprecated function


private:
  char *command;               // ptr to current command
  string line,copy;            // input line string
  int maxline, maxcopy;        // max lengths of char strings
  int maxarg;                  // max # of args in arg

  int numtriple(char *);                     // count number of triple quotes
  double precedence(char);                   // find precedence of operators.
  Var applyOp(Var, Var, char);              // perform arithmetic operations.
  bool is_operator(char);                    // check if is an operator
  bool is_math_char(char);                   // check if the character is either of +-/*()
  Var evaluate_function(string , string); // evaluate function with argument
  string remove_whitespace(string);          // remove white spaces from string
  int dimension(vector<string>);             // set the dimension of the simulation domain
  int region(vector<string>);
  int solid(vector<string>);
  int EOS(vector<string>);
  int dump(vector<string>);
  int group_command(vector<string>);
  int log(vector<string>);
  int method_modify(vector<string>);
  int fix(vector<string>);

 public:
  typedef void (*CommandCreator)(MPM *,vector<string>);
  typedef map<string,CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;

 protected:
  template <typename T> static void command_creator(MPM *,vector<string>);

};

#endif
