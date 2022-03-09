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

#include <pointers.h>
#include <var.h>
#include <vector>
#include <map>


/*! This class reads the input file and execute the instructions present in it line by line.
 * 
 * This class is the back bone of the code as it is the interface between the user and the code (through the input file).
 * Its member functions are design in order to be able to generate variables and execute functions present in the input file.
 * This class has Pointers as parent in order to be accessible from anywhere within the MPM class.
 */
class Input : protected Pointers {
  friend class Error;
public:
  int narg;                    ///< Number of command args
  char **arg;                  ///< Parsed args for command

  map<string, class Var> *vars; ///< List of global user variables.
                                ///< These variables are those created by the user through the input file.

  // functions
  Input(class MPM *, int, char **);
  ~Input();
  void file();                 ///< Reads the input file line by line and pass it to parsev().
  class Var parsev(string);    ///< Parse an input text line.
  double parse(string);        ///< Deprecated function


private:
  int me;                      ///< Proc ID
  char *command;               ///< Pointer to current command
  string line, copy;           ///< input line string
  int line_number;             ///< line number of the input file being processed
  int maxline, maxcopy;        ///< max lengths of char strings
  int maxarg;                  ///< max number of args in arg

  int numtriple(char *);                     ///< Counts the number of triple quotes
  double precedence(const string);           ///< Finds precedence of operators.
  Var applyOp(Var, const string, Var);       ///< Performs arithmetic operations.
  bool is_operator(char);                    ///< Checks if op is an operator. Return true or false.
  bool is_math_char(char);                   ///< Checks if the character is either of +-/*()
  Var evaluate_function(string , string);    ///< Evaluates the user function with argument
  string remove_whitespace(string);          ///< Removes white spaces from string
  int dimension(vector<string>);             ///< Sets the dimension of the simulation domain
  int axisymmetric(vector<string>);          ///< Sets the simulation to axisymmetric
  int region(vector<string>);                ///< Creates a region
  int solid(vector<string>);                 ///< Creates a solid
  int add_EOS(vector<string>);               ///< Adds an Equation of State (EOS)
  int add_strength(vector<string>);          ///< Adds a strength
  int add_material(vector<string>);          ///< Adds a material
  int add_damage(vector<string>);            ///< Adds a damage law
  int add_temperature(vector<string>);       ///< Adds a temperature law
  int dump(vector<string>);                  ///< Dump a snapshot to be seen with Ovito
  int group_command(vector<string>);         ///< Creates a group of particles and/or nodes
  int set_output(vector<string>);            ///< Sets the frequency of console and logfile outputs
  int log_modify(vector<string>);            ///< Modifies the varibles to be displayed on the console and saved in the log file.
  int method(vector<string>);                ///< Specifies the MPM method used (tlmpm, ulmpm, ...)
  int scheme(vector<string>);                ///< Creates a scheme, USL, or MUSL.
  int fix(vector<string>);                   ///< Creates a fix.
  int delete_fix(vector<string>);            ///< Deletes a fix.
  int compute(vector<string>);               ///< Creates a compute.
  int delete_compute(vector<string>);        ///< Deletes a compute.
  int set_dt_factor(vector<string>);         ///< Sets the factor to be applied to the CFL timestep
  int set_dt(vector<string>);                ///< Sets the timestep
  class Var value(vector<string>);           ///< Returns the current value of a user variable.
  int plot(vector<string>);                  ///< Add a curve to be plotted.
  int save_plot(vector<string>);             ///< Save the plot as ...
  int print(vector<string>);                 ///< Print the details of a user-variable.
  int create_domain(vector<string>);         ///< Deprecated.
  int restart(vector<string>);               ///< Create a restart policy

  bool protected_variable(string);           ///< Checks if the variables is protected or not.

  vector<string> protected_vars;             ///< List of protected variables.

 public:
  typedef class Var (*CommandCreator)(MPM *,vector<string>);
  typedef map<string,CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;           ///< Map of the commands listed in style_command.h

 protected:
  template <typename T> static class Var command_creator(MPM *,vector<string>);

};

#endif
