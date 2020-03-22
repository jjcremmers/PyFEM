lc = 5e-2;

Point(1) = { 0,  0, 0, lc} ;
Point(2) = {.1,  0, 0, 0.2*lc} ;
Point(3) = {.5,  0, 0, lc} ;
Point(4) = {.5, 1., 0, lc} ;
Point(5) = { 0, 1., 0, lc} ;
Point(6) = { 0, .1, 0, 0.2*lc} ;

Line(1) = {2,3} ;
Line(2) = {3,4} ;
Line(3) = {4,5} ;
Line(4) = {5,6} ;

Circle(5) = {6,1,2};


Line Loop(1) = {1,2,3,4,5};

Plane Surface(1) = {1};

Physical Surface("My surface") = {1};

