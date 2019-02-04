#include "mpm.h"
#include "input.h"
#include "domain.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stack>

#define DELTALINE 256
#define DELTA 4

using namespace std;


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
  bool ignore = false;

  istream is(infile);
  while (is) {
    char c = char(is.get());
    if (c != '\n') {
      if (c == '#') ignore = true; // ignore everything after #

      if (!ignore)
	line.append(&c,1);
    } else {
      ignore = false;
      cout << line << endl;
      cout << "gives: " << parse(line) << endl;
      line.clear();
    }
  }  
}

// Function to find precedence of  
// operators. 
double Input::precedence(char op){ 
    if(op == '+'||op == '-') return 1;
    if(op == '*'||op == '/') return 2;
    if(op == '^') return 3;
    return 0;
} 
  
// Function to perform arithmetic operations. 
double Input::applyOp(double a, double b, char op){ 
    switch(op){ 
    case '+': return a + b;
    case '-': return a - b;
    case '*': return a * b;
    case '/': return a / b;
    case '^': return pow(a,b);
    case '(':
      printf("Error: unmatched parenthesis (\n");
      exit(1);
    default:
      printf("Error: unknown operator %c\n", op);
      exit(1);
    } 
}


bool Input::is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  return false;
}

// check if op is either of +-/*()
bool Input::is_math_char(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  if (op=='^') return true;
  if (op=='(') return true;
  if (op==')') return true;
  if (op=='=') return true;
  return false;
}

// evaluate function func with argument arg:
double Input::evaluate_function(string func, string arg){
  cout << "Evaluate function " << func << " with argument: " << arg << endl;


  if (func.compare("dimension") == 0) return (double) dimension(arg);
  if (func.compare("region") == 0) return (double) region(arg);
  else if (func.compare("exp") == 0) return (double) exp(parse(arg));
  cout << "Error: Unknown function " << func << endl;
  exit(1);
}

// remove white spaces from string
string Input::remove_whitespace(string str){
  string str_;

  for(int i=0; i<str.length(); i++){
    if( str[i] != ' ') str_.append(&str[i],1); // Add the non-whitespace character to str_
  }
  return str_;
}


double Input::parse(string str)
{
  // stack to store integer values. 
  stack <double> values; 
      
  // stack to store operators. 
  stack <char> ops;

  // stack to store functions.
  stack <string> funcs;

  string returnvar;

  str = remove_whitespace(str);
  cout << "New string: " << str << endl;

  bool negative = false; // flag that indicates that the next value will have to be multiplied by -1

  for (int i=0; i<str.length(); i++) {
    if (isdigit(str[i])) {

      string number;
      
      // The number could be longer than one digit:
      int j;
      for (j=1; j<str.length()-i;j++) {
	if (!isdigit(str[i+j]) && str[i+j]!='.') break;
      }

      number.append(&str[i],j);
      i += j-1; 

      if (negative) {
      	values.push((-1)*stof(number));
      	negative = false;
      } else {
      	values.push(stof(number));
      }
      
      cout << "Pushed number: " << values.top() << " to values stack" << endl;
    }

    // If the first character in the righ-hand-side is (, push it to ops
    else if (str[i] == '(') {
      cout << "Pushed \'(\' to ops stack" << endl;
      ops.push(str[i]);

      if (i+1 < str.length() && str[i+1] == '-') {
	negative = true;
	i++;
      }
    }

    // Closing brace encountered, solve
    // entire brace.
    else if (str[i] == ')') {
      cout << "Found \')\'" << endl;

      if (ops.empty() || values.size() < 2) {
	printf("Error, unmatched parenthesis )\n");
	exit(1);
      }

      while(ops.top() != '(')
	{

	  double val2 = values.top();
	  values.pop();

	  double val1 = values.top();
	  values.pop();

	  char op = ops.top();
	  ops.pop();
	  printf("%f %f %c\n", val1, val2, op);

	  values.push(applyOp(val1, val2, op));

	  if (ops.empty()) {
	    printf("Error, unmatched parenthesis )\n");
	    exit(1);
	  }
	}

      // pop opening brace. 
      ops.pop();
    }


    else if (is_operator(str[i])){
      char new_op = str[i];
      printf("found operator %c\n", new_op);

      if (values.empty()) {
	if (new_op == '-') {
	  negative = true;
	  continue;
	}
      }

      if (i+1 >= str.length()) {
	printf("Error: end-of-line character detected after operator %c\n", new_op);
	exit(1);
      }

      if (i+1 < str.length() && str[i+1] == '\0') {
	printf("Error: end-of-line character detected after operator %c\n", new_op);
	exit(1);
      }

      else if (i+1 < str.length() && str[i+1] == '*') {
	new_op = '^';
	i++;
      }

      else if (i+1 < str.length() && is_operator(str[i+1]) && str[i+1] != '-') {
	printf("Error: unknown operator sequence %c%c\n", new_op, str[i+1]);
	exit(1);
      }


      // While top of 'ops' has same or greater
      // precedence to current token, which
      // is an operator. Apply operator on top
      // of 'ops' to top two elements in values stack.
      while(!ops.empty() &&
	    precedence(ops.top()) >= precedence(new_op)){
	double val2 = values.top();
	values.pop();

	double val1 = values.top();
	values.pop();

	char op = ops.top();
	ops.pop();
	printf("%f %f %c\n", val1, val2, op);

	values.push(applyOp(val1, val2, op));
      } 
              
      // Push current token to 'ops'. 
      ops.push(new_op);

      if (i+1 < str.length() && str[i+1] == '-') {
	negative = true;
	i++;
      }
    }


    else {
      string word;
      
      // Check how long is the word
      int j;
      for (j=1; j<str.length()-i;j++) {
	if (is_math_char(str[i+j])) break;
      }

      word.append(&str[i],j);
      i += j-1; 

      cout << "Found keyword: " << word << endl;

      // Check if word is a variable:
      map<string, double>::iterator it;
      it = variables->find(word);

      if (it != variables->end()){
	// word is a variable
	cout << "word is a variable\n";
	if (i+1 < str.length() && str[i+1] == '=') {

	  if (!values.empty() || !ops.empty() ) {
	    printf("Error: I do not understand when '=' is located in the middle of an expression\n");
	    exit(1);
	  }

	  else {
	    returnvar = word;
	    cout << "The computed value will be stored in " <<  returnvar << endl;
	    i++;
	  }
	}

	else if (i+1 >= str.length() && values.empty() && ops.empty() ) {
	  if (negative) {
	    return (-1)*(*variables)[word];
	  }
	  else {
	    return (*variables)[word];
	  }
	}
	
	else {
	  if (negative) {
	    values.push((-1)*(*variables)[word]);
	    negative = false;
	  } else {
	    values.push((*variables)[word]);
	  }
	  cout << "push " << word << "=" << values.top() << " to values\n";
	}
      }

      else if (i+1 < str.length() && str[i+1]=='=' && values.empty() && ops.empty()){
	// Check if there is an '=':
	cout << "Check if there is a =\n"; 
	returnvar = word;
	cout << "The computed value will be stored in " <<  returnvar << endl;
	i++;
      }

      else if (i+1 < str.length() && str[i+1]=='(') {
	// It's a function:
	i+=2;

	// Extract the argument:
	int k = 0;
	while(str[i+k]!=')'){
	  k++;
	  if (i+k > str.length()) {
	    cout << "Error: Unbalanced parenthesis '('" << endl;
	    exit(1);
	  }
	}
	string arg;
	arg.append(&str[i],k);
	cout << "Found function " << word << " with argument: " << arg << endl;
	i += k;
	values.push(evaluate_function(word, arg));
      }

      else {
	cout << word << " is unknown\n";
	exit(1);
      }
    }
  }

  // Entire expression has been parsed at this
  // point, apply remaining ops to remaining
  // values.
  while(!ops.empty()){
    double val2 = values.top();
    values.pop();

    double val1 = values.top();
    values.pop();

    char op = ops.top();
    ops.pop();
                  
    values.push(applyOp(val1, val2, op)); 
  }

  // Top of 'values' contains result, return it.
  if (values.empty()) return -1;
  else {
    printf("value=%f\n", values.top());
    if (!returnvar.empty()) {
      (*variables)[returnvar] = values.top();
    }
    return values.top();
  }
}


int Input::dimension(string arg){

  int dim = (int) parse(arg);


  if (dim != 2 && dim != 3) {
    cout << "Error: dimension argument: " << dim << endl;
    exit(1);
  }
  else domain->dimension = dim;

  cout << "Set dimension to " << dim << endl;
  return 0;
}


int Input::region(string str){
  // Separate arguments:
  
  int j = 0;
  int start = 0;
  for (int i=0; i<str.length(); i++) {
    // Locate comas.
    if (str[i] == ',' || i==str.length()-1)  {
      if (i==start) {
	cout << "Error: missing argument" << endl;
	exit(1);
      }

      args.resize(args.size()+1);
      args.back().append(&str[start], i - start + (i==str.length()-1) );
      cout << "Received argument " << j+1 << " :" << args.back() << endl;
      start = i+1;
      j++;
    }
  }

  domain->add_region(args);
  return 0;
}
