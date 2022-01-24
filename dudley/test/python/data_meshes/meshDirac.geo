Mesh.MshFileVersion=4.4;
// basic box with 5 electrodes
WX =  2.00000e+02;
WY =  1.00000e+02;
WZ = -1.00000e+02;
wx = WX/2.;
wy = WY/2.;

// element sizes
mshC = 1.600000e+01;
mshE = 1.000000e+00;

// Box top points
Point(1) = { -wx, -wy, 0.0, mshC};
Point(2) = {  wx, -wy, 0.0, mshC};
Point(3) = { -wx,  wy, 0.0, mshC};
Point(4) = {  wx,  wy, 0.0, mshC};

// Big Box bottom points
Point(5) = { -wx, -wy, WZ, mshC};
Point(6) = {  wx, -wy, WZ, mshC};
Point(7) = { -wx,  wy, WZ, mshC};
Point(8) = {  wx,  wy, WZ, mshC};

// all lines in positive axes directions
// Box top
Line(1) = {1,2} ;
Line(2) = {3,4} ;
Line(3) = {1,3} ;
Line(4) = {2,4} ;
// Box bottom
Line(5) = {5,6} ;
Line(6) = {7,8} ;
Line(7) = {5,7} ;
Line(8) = {6,8} ;
// Box sides
Line(9) = {1,5} ;
Line(10) = {2,6} ;
Line(11) = {3,7} ;
Line(12) = {4,8} ;

// Big Box surfaces
// outward normals!!!
Line Loop(25) = {1,4,-2,-3} ;   // top    z = 0
Line Loop(26) = {7,6,-8,-5} ;   // bottom z = WZ
Line Loop(27) = {9,5,-10,-1} ;  // front  y = -wy
Line Loop(28) = {2,12,-6,-11} ; // back   y =  wy
Line Loop(29) = {3,11,-7,-9} ;  // left   x = -wx
Line Loop(30) = {10,8,-12,-4} ; // right  x =  wx

Plane Surface(1) = {25};  // top    z = 0
Plane Surface(2) = {26};  // bottom z = WZ
Plane Surface(3) = {27};  // front  y = -wy
Plane Surface(4) = {28};  // back   y =  wy
Plane Surface(5) = {29};  // left   x = -wx
Plane Surface(6) = {30};  // right  x =  wx

k = newp;
// electrodes  
Point(k) = { -6.400000e+01, 0.000000e+00, 0.0, mshE};
Point{k} In Surface{1};
Point(k+1) = { -3.200000e+01, 0.000000e+00, 0.0, mshE};
Point{k+1} In Surface{1};
Point(k+2) = { 0.000000e+00, 0.000000e+00, 0.0, mshE};
Point{k+2} In Surface{1};
Point(k+3) = { 3.200000e+01, 0.000000e+00, 0.0, mshE};
Point{k+3} In Surface{1};
Point(k+4) = { 6.400000e+01, 0.000000e+00, 0.0, mshE};
Point{k+4} In Surface{1};
// Surface loops and volumes

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Volume("volume")  = { 1 } ;
Physical Surface("top")    = { 1 };
Physical Surface("bottom") = { 2 };
Physical Surface("front")  = { 3 };
Physical Surface("back")   = { 4 };
Physical Surface("left")   = { 5 };
Physical Surface("right")  = { 6 };

Physical Point("e0") = {k};
Physical Point("e1") = {k+1};
Physical Point("e2") = {k+2};
Physical Point("e3") = {k+3};
Physical Point("e4") = {k+4};





