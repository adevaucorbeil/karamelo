#include "mpm.h"
#include "input.h"
#include <iostream>
#include <fstream>
#include <string.h>

#define DELTALINE 256
#define DELTA 4

Input::Input(MPM *mpm, int argc, char **argv) : Pointers(mpm)
{
  maxline = maxcopy = 0;
  maxarg = 0;
  arg = NULL;
}


Input::~Input()
{
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int m,n;

  while (1) {

    // read a line from input script
    // n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line

    m = 0;
    while (1) {
      if (maxline-m < 2) reallocate(line,maxline,0);

      // end of file reached, so break
      // n == 0 if nothing read, else n = line with str terminator

      if (fgets(&line[m],maxline-m,infile) == NULL) {
	if (m) n = strlen(line) + 1;
	else n = 0;
	break;
      }

      // continue if last char read was not a newline
      // could happen if line is very long

      m = strlen(line);
      if (line[m-1] != '\n') continue;

      // continue reading if final printable char is & char
      // or if odd number of triple quotes
      // else break with n = line with str terminator

      m--;
      while (m >= 0 && isspace(line[m])) m--;
      if (m < 0 || line[m] != '&') {
	if (numtriple(line) % 2) {
	  m += 2;
	  continue;
	}
	line[m+1] = '\0';
	n = m+2;
	break;
      }
    }
    if (n==0) break;
    parse();
    if (command == NULL) continue;
  }
}


/* ----------------------------------------------------------------------
   rellocate a string
   if n > 0: set max >= n in increments of DELTALINE
   if n = 0: just increment max by DELTALINE
------------------------------------------------------------------------- */

void Input::reallocate(char *&str, int &max, int n)
{
  if (n) {
    while (n > max) max += DELTALINE;
  } else max += DELTALINE;

  str = (char *) realloc(str,max*sizeof(char));
}

/* ----------------------------------------------------------------------
   return number of triple quotes in line
------------------------------------------------------------------------- */

int Input::numtriple(char *line)
{
  int count = 0;
  char *ptr = line;
  while ((ptr = strstr(ptr,"\"\"\""))) {
    ptr += 3;
    count++;
  }
  return count;
}


void Input::parse()
{
  // duplicate line into copy string to break into words

  int n = strlen(line) + 1;
  if (n > maxcopy) reallocate(copy,maxcopy,n);
  strcpy(copy,line);

  int id_last_space = -1;
  int id_space = -1;
  int nword = 0;

  for (int i=0; i<n; i++){
    if (copy[i] == '#') {
      copy[i] = '\0';
      break;
    }
  }

  char *next;
  command = nextword(copy,&next);
  printf("command = %s\n", command);
  if (command == NULL) return;

  narg = 0;
  char *ptr;
  ptr = next;
  while (ptr) {
    if (narg == maxarg) {
      maxarg += DELTA;
      //printf("realloc\n");
      arg = (char **) realloc(arg,maxarg*sizeof(char *));
      if (arg == NULL) {
	printf("arg == NULL\n");
	exit(1);
      }
    }
    arg[narg] = nextword(ptr,&next);
    if (!arg[narg]) break;
    narg++;
    ptr = next;
  }
  for (int i=0; i<narg; i++) {
    printf("%s\n",arg[i]);
  }
  printf("\n");
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   treat text between single/double/triple quotes as one arg
   matching quote must be followed by whitespace char if not end of string
   strip quotes from returned word
   return ptr to start of word or NULL if no word in string
   also return next = ptr after word
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  // start = first non-whitespace char

  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return NULL;

  // if start is single/double/triple quote:
  //   start = first char beyond quote
  //   stop = first char of matching quote
  //   next = first char beyond matching quote
  //   next must be NULL or whitespace
  // if start is not single/double/triple quote:
  //   stop = first whitespace char after start
  //   next = char after stop, or stop itself if stop is NULL

  if (strstr(start,"\"\"\"") == start) {
    stop = strstr(&start[3],"\"\"\"");
    // if (!stop) //error->all(FLERR,"Unbalanced quotes in input line");
    start += 3;
    *next = stop+3;
    // if (**next && !isspace(**next))
    //   error->all(FLERR,"Input line quote not followed by whitespace");
  } else if (*start == '"' || *start == '\'') {
    stop = strchr(&start[1],*start);
    // if (!stop) error->all(FLERR,"Unbalanced quotes in input line");
    start++;
    *next = stop+1;
    // if (**next && !isspace(**next))
    //   error->all(FLERR,"Input line quote not followed by whitespace");
  } else {
    stop = &start[strcspn(start," \t\n\v\f\r")];
    if (*stop == '\0') *next = stop;
    else *next = stop+1;
  }

  // set stop to NULL to terminate word

  *stop = '\0';
  return start;
}
