/* -*- c++ -*- ----------------------------------------------------------
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

#ifndef MPM_INPUT_H
#define MPM_INPUT_H

#include "pointers.h"
#include "var.h"
#include <vector>
#include <map>


class Input : protected Pointers {
  friend class Error;
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
  int me;                      // proc ID
  char *command;               // ptr to current command
  string line,copy;            // input line string
  int line_number;             // line number of the input file being processed
  int maxline, maxcopy;        // max lengths of char strings
  int maxarg;                  // max # of args in arg

  int numtriple(char *);                     // count number of triple quotes
  double precedence(string);                   // find precedence of operators.
  Var applyOp(Var, Var, string);              // perform arithmetic operations.
  bool is_operator(char);                    // check if is an operator
  bool is_math_char(char);                   // check if the character is either of +-/*()
  Var evaluate_function(string , string); // evaluate function with argument
  string remove_whitespace(string);          // remove white spaces from string
  int dimension(vector<string>);             // set the dimension of the simulation domain
  int axisymmetric(vector<string>);          // set the simulation to axisymmetric
  int region(vector<string>);
  int solid(vector<string>);
  int add_EOS(vector<string>);
  int add_strength(vector<string>);
  int add_material(vector<string>);
  int add_damage(vector<string>);
  int add_temperature(vector<string>);
  int dump(vector<string>);
  int group_command(vector<string>);
  int set_output(vector<string>);
  int log_modify(vector<string>);
  int method(vector<string>);
  int scheme(vector<string>);
  int fix(vector<string>);
  int set_dt_factor(vector<string>);
  int set_dt(vector<string>);
  class Var value(vector<string>);
  int plot(vector<string>);
  int save_plot(vector<string>);
  int print(vector<string>);
  int create_domain(vector<string>);

  bool protected_variable(string);

  vector<string> protected_vars;

 public:
  typedef class Var (*CommandCreator)(MPM *,vector<string>);
  typedef map<string,CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;

 protected:
  template <typename T> static class Var command_creator(MPM *,vector<string>);

};

#endif
