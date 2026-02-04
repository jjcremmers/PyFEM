lc = 0.5e-3;

Point(1)  = { 0.0, 0.0, 0.0, lc} ;
Point(2)  = { 2.0e-2, 0.0, 0.0, lc} ;
Point(3)  = { 2.0e-2, 2.0e-2, 0.0, lc} ;
Point(4)  = { 0.0, 2.0e-2, 0.0, lc} ;

Point(5)  = { 1.5e-3, 6.e-3, 0.0, lc} ;
Point(6)  = { 6.5e-3, 6.e-3, 0.0, lc} ;
Point(7)  = { 11.5e-3, 6.e-3, 0.0, lc} ;

Point(8)  = { 9.0e-3, 13.5e-3, 0.0, lc} ;
Point(9)  = { 14.0e-3, 13.5e-3, 0.0, lc} ;
Point(10) = { 19.0e-3, 13.5e-3, 0.0, lc} ;


Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Circle(5) = {5,6,7};
Circle(6) = {7,6,5};

Circle(7) = {8,9,10};
Circle(8) = {10,9,8};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6};
Line Loop(3) = {7,8};

Plane Surface(1) = {1,-2,-3};
Plane Surface(2) = {2};
Plane Surface(3) = {3};


Physical Surface("Epoxy")  = {1};


Physical Curve("Left")      = (4);
Physical Curve("Right")     = (2);
Physical Curve("Top")       = (3);
Physical Curve("Bottom")    = (1);
Physical Point("SinglePoint") = (3);
