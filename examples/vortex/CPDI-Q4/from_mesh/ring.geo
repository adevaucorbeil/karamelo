ra=0.75;
rb=1.25;
h=1;

n1 = 10; 
n2 = 7;

Point(1)={0,0,0,h};

Point(2)={ra,0,0,h};
Point(3)={0,ra,0,h};
Point(4)={-ra,0,0,h};
Point(5)={0,-ra,0,h};

Point(22)={rb,0,0,h};
Point(32)={0,rb,0,h};
Point(42)={-rb,0,0,h};
Point(52)={0,-rb,0,h};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};

Circle(12)={22,1,32};
Circle(22)={32,1,42};
Circle(32)={42,1,52};
Circle(42)={52,1,22};

Line(5) = {2,22};
Line(6) = {3,32};
Line(7) = {4,42};
Line(8) = {5,52};

Line Loop(1) = {5,12,-6,-1};
Line Loop(2) = {6,22,-7,-2};
Line Loop(3) = {7,32,-8,-3};
Line Loop(4) = {8,42,-5,-4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Transfinite Line{1,2,3,4,12,22,32,42}=n1;
Transfinite Line{5,6,7,8}=n2;

Transfinite Surface{1} = {2,22,32,3};
Transfinite Surface{2} = {3,32,42,4};
Transfinite Surface{3} = {4,42,52,5};
Transfinite Surface{4} = {5,52,22,2};

Recombine Surface{1,2,3,4};

Physical Surface(888) = {1,2,3,4};






