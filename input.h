/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_INPUT_H
#define MPM_INPUT_H

#include "pointers.h"

class Input : protected Pointers {
public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command

  // functions
  Input(class MPM *, int, char **);
  ~Input();
  void file();                 // process all input


private:
  char *command;               // ptr to current command
  string line,copy;            // input line string
  int maxline, maxcopy;        // max lengths of char strings
  int maxarg;                  // max # of args in arg

  void reallocate(char *&, int &, int);  // reallocate a char string
  int numtriple(char *);                 // count number of triple quotes
  void parse();                          // parse an input text line
  char *nextword(char *, char **);       // find next word in string with quotes
};

#endif
