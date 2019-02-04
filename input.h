/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_INPUT_H
#define MPM_INPUT_H

#include "pointers.h"
#include <vector>

class Input : protected Pointers {
public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  vector<string> args;

  // functions
  Input(class MPM *, int, char **);
  ~Input();
  void file();                 // process all input
  double parse(string);        // parse an input text line


private:
  char *command;               // ptr to current command
  string line,copy;            // input line string
  int maxline, maxcopy;        // max lengths of char strings
  int maxarg;                  // max # of args in arg

  int numtriple(char *);                     // count number of triple quotes
  double precedence(char);                   // find precedence of operators.
  double applyOp(double , double , char );   // perform arithmetic operations.
  bool is_operator(char);                    // check if is an operator
  bool is_math_char(char);                   // check if the character is either of +-/*()
  double evaluate_function(string , string); // evaluate function with argument
  string remove_whitespace(string);          // remove white spaces from string
  int dimension(string);                     // set the dimension of the simulation domain
  int region(string);

};

#endif
