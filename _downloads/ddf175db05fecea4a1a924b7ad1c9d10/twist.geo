a  = 1.0;
b  = 3.0;
lc = 1.0;

Point(1)  = { -a , -b , 0.0 , lc };
Point(2)  = {  a , -b , 0.0 , lc };
Point(3)  = {  a , -a , 0.0 , lc };
Point(4)  = {  b , -a , 0.0 , lc };
Point(5)  = {  b ,  a , 0.0 , lc };
Point(6)  = {  a ,  a , 0.0 , lc };
Point(7)  = {  a ,  b , 0.0 , lc };
Point(8)  = { -a ,  b , 0.0 , lc };
Point(9)  = { -a ,  a , 0.0 , lc };
Point(10) = { -b ,  a , 0.0 , lc };
Point(11) = { -b , -a , 0.0 , lc };
Point(12) = { -a , -a , 0.0 , lc };

Line(1)  = {1,2} ;
Line(2)  = {2,3} ;
Line(3)  = {3,4} ;
Line(4)  = {4,5} ;
Line(5)  = {5,6} ;
Line(6)  = {6,7} ;
Line(7)  = {7,8} ;
Line(8)  = {8,9} ;
Line(9)  = {9,10} ;
Line(10) = {10,11} ;
Line(11) = {11,12} ;
Line(12) = {12,1} ;

Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};

Plane Surface(1) = {1};

e() = Extrude{{0,0,10},{0,0,1},{0,0,0},Pi/2} { Surface{1};Layers{10};Recombine;
};

Physical Surface("Bottom") = {1};
Physical Surface("Top")    = {74};
Physical Volume("Body")    = {1};

