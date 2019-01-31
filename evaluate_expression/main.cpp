#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <bits/stdc++.h> 

using namespace std; 

map<string, int> variables;

// Function to find precedence of  
// operators. 
int precedence(char op){ 
    if(op == '+'||op == '-') 
    return 1; 
    if(op == '*'||op == '/') 
    return 2; 
    return 0; 
} 
  
// Function to perform arithmetic operations. 
int applyOp(int a, int b, char op){ 
    switch(op){ 
    case '+': return a + b; 
    case '-': return a - b; 
    case '*': return a * b; 
    case '/': return a / b;
    case '(':
      printf("Error: unmatched parenthesis (\n");
      exit(1);
    default:
      printf("Error: unknown operator %c\n", op);
      exit(1);
    } 
}

void extract_lhs(char *str, char *&lhs, char *&rhs)
{
  size_t size;

  size =  strcspn(str,"+-*/()=");
  //if (size == strlen(str)) return double()
  lhs = NULL;
  if (size > 0) lhs = (char*) malloc(size+1);
  rhs = strpbrk(str,"+-*/()=");
  if (rhs != NULL) {
    //rhs++;
    rhs = &rhs[strspn(rhs," ")];
  }
  if (size > 0) {
    strncpy(lhs,str,size);
    lhs[(int) size] = '\0'; // Manually add null character at the end of the string
  }
  printf("extract_lhs: str=%s, lhs=%s, rhs=%s\n", str, lhs, rhs);
}

bool is_operator(char op){
  if (op=='+') return true;
  if (op=='-') return true;
  if (op=='*') return true;
  if (op=='/') return true;
  return false;
}




int evaluate(char *str)
{
  // stack to store integer values. 
  stack <int> values; 
      
  // stack to store operators. 
  stack <char> ops; 

  string returnvar;
  
  char *e, *lhs, *rhs;
  e = NULL;
  lhs = NULL;
  rhs = NULL;
  rhs = (char*) malloc(strlen(str)+1);
  strcpy(rhs, str);

  //for (int i=0;i<3;i++) {
  while (rhs!=NULL){
    if(isdigit(rhs[0])){
      e = rhs;
      rhs = NULL;
      extract_lhs(e, lhs, rhs);
      printf("push %d to values\n", atoi(lhs));
      values.push(atoi(lhs));
    }

    // If the first character in the righ-hand-side is (, push it to ops
    else if(rhs[0] == '('){
      printf("found (\n");
      ops.push(rhs[0]);
      rhs++;
      rhs = &rhs[strspn(rhs," ")]; // Remove whitespace
      printf("( lhs=%s, rhs=%s\n", lhs, rhs);
    }

    // Closing brace encountered, solve  
    // entire brace. 
    else if(rhs[0] == ')') {
      printf("found )\n");
      while(ops.top() != '(') 
	{
	  int val2 = values.top(); 
	  values.pop(); 
          
	  int val1 = values.top(); 
	  values.pop(); 
                  
	  char op = ops.top(); 
	  ops.pop(); 
	  printf("%d %d %c\n", val1, val2, op);
                  
	  if (ops.empty()) {
	    printf("Error, unmatched parenthesis )\n");
	    exit(1);
	  }
	  values.push(applyOp(val1, val2, op)); 
	} 
              
      // pop opening brace. 
      ops.pop();
      rhs++;
      if (rhs[0] == '\0') rhs = NULL;
      else rhs = &rhs[strspn(rhs," ")]; // Remove whitespace
      printf(") lhs=%s, rhs=%s\n", lhs, rhs);
    } 

    else if (is_operator(rhs[0])){
      printf("found operator %c\n", rhs[0]);
      // While top of 'ops' has same or greater  
      // precedence to current token, which 
      // is an operator. Apply operator on top  
      // of 'ops' to top two elements in values stack. 
      while(!ops.empty() && precedence(ops.top()) 
	    >= precedence(rhs[0])){ 
	int val2 = values.top(); 
	values.pop(); 
                  
	int val1 = values.top(); 
	values.pop(); 
        
	char op = ops.top(); 
	ops.pop();
	printf("%d %d %c\n", val1, val2, op);
                  
	values.push(applyOp(val1, val2, op)); 
      } 
              
      // Push current token to 'ops'. 
      ops.push(rhs[0]);
      rhs++;
      rhs = &rhs[strspn(rhs," ")]; // Remove whitespace
    }

    else {
      e = rhs;
      rhs = NULL;
      extract_lhs(e, lhs, rhs);

      
      printf("lhs=%s, rhs=%s\n", lhs, rhs);

      // Check if lhs is a variable:
      map<string, int>::iterator it;
      it = variables.find(lhs);
      if (it != variables.end()){
	values.push(variables[(string) lhs]);
	printf("push %s=%d to values\n", lhs, values.top());
      } else {
	if (rhs[0]=='=' && values.empty() && ops.empty()){
	  returnvar = string(lhs);
	  rhs++;
	  rhs = &rhs[strspn(rhs," ")]; // Remove whitespace
	  cout << "The computed value will be stored in " <<  returnvar << endl;
	} else {
	  printf("%s is unknown\n", lhs);
	  exit(1);
	}
      }
    }
  }

  // Entire expression has been parsed at this 
  // point, apply remaining ops to remaining 
  // values. 
  while(!ops.empty()){ 
    int val2 = values.top(); 
    values.pop(); 
                  
    int val1 = values.top(); 
    values.pop(); 
                  
    char op = ops.top(); 
    ops.pop(); 
    printf("%d %d %c\n", val1, val2, op);
                  
    values.push(applyOp(val1, val2, op)); 
  } 
      
  // Top of 'values' contains result, return it.
  if (!returnvar.empty()) {
    variables[returnvar] = values.top();
    cout << returnvar << "=" <<  variables[returnvar] << endl;
  }
  return values.top(); 
}






int main(){
  variables["ab"] = 10;
  
  double result;
  char *expr = "x= 2 + (11*ab) + 10";

  printf("%s = %d\n", expr, evaluate(expr));
}


