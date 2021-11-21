#pragma once

#include <Mat2.h>

void svd(const Mat2 &A,
         Mat2 &U,
         Mat2 &sig,
         Mat2 &V) {
  double x = A(0, 0) + A(1, 1), y = A(1, 0) - A(0, 1);
  double scale = 1/sqrt(x*x + y*y);
  double c = x*scale;
  double s = y*scale;
  Mat2 R(c, -s, s, c);
  Mat2 S = R.transpose()*A;
  double s1 = 0, s2 = 0;
  if (abs(S(0, 1)) < 1e-5) {
    c = 1;
    s = 0;
  }
  else {
   double tao = 0.5*(S(0, 0) - S(1, 1));
   double w = sqrt(tao*tao + S(0, 1)* S(0, 1));
   double t = 0;
   if (tao > 0)
     t = S(0, 1)/(tao + w);
   else
     t = S(0, 1)/(tao - w);
   c = 1/sqrt(t*t + 1);
   s = -t*c;
   s1 = c*c*S(0, 0) - 2*c*s*S(0, 1) + s*s*S(1, 1);
   s2 = s*s*S(0, 0) + 2*c*s*S(0, 1) + c*c*S(1, 1);
  }
       
  if (s1 < s2) {
    double tmp = s1;
    s1 = s2;
    s2 = tmp;
    V = Mat2(-s, c, -c, -s);
  }
  else
    V = Mat2(c, s, -s, c);
  U = R*V;
  sig = Mat2(s1, 0, 0, s2);
}
