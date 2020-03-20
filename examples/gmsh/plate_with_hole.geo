lc = 1e-2;

Point(1) = { 0,  0, 0, lc} ;
Point(2) = {.2,  0, 0, lc} ;
Point(3) = {.2, .2, 0, lc} ;
Point(4) = { 0, .2, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Point(5) = { .07 , 0.1  , 0 , 0.5*lc };
Point(6) = { .1  , 0.13 , 0 , 0.5*lc };
Point(7) = { .13 , 0.1  , 0 , 0.5*lc };
Point(8) = { .10 , 0.07 , 0 , 0.5*lc };
Point(9) = { .1  , 0.1  , 0 , 0.5*lc };

Circle(5) = {5,9,6};
Circle(6) = {6,9,7};
Circle(7) = {7,9,8};
Circle(8) = {8,9,5};


Line Loop(1) = {4,1,-2,3} ;
Line Loop(2) = {5,6,7,8};

Plane Surface(1) = {1};

Physical Surface("My surface") = {1};

