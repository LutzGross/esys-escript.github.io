// Example of the use of unstructured mesh in esys.finley
//
// Run this to generate mesh:
//
//      gmsh -2 -format msh22 -o inclusion.msh inclusion.geo
//
Mesh.MshFileVersion = 2.2; // Version of the MSH file format to use
//The geometry is a rectangle with an inclusion.
//
// ... Dimensions:
L0 = 200.;
L1 = 100.;
X0_c = 80;
X1_c = 60;
R_c = 20.;ignosre

// Mesh sizes:
SizeInter = 2.;
SizeOuter = 10.;

// ... Points: Outer Boundary
Point(1) = {0., 0., 0., SizeOuter};
Point(2) = {L0, 0., 0., SizeOuter};
Point(3) = {0., L1, 0., SizeOuter};
Point(4) = {L0, L1, 0., SizeOuter};
//... Points: Inner Boundary
// ..... Center of Circle:
Point(5) = {X0_c , X1_c, 0., SizeInter};
// ..... Points on the perimeter:
Point(6) = {X0_c, X1_c - R_c, 0., SizeInter};
Point(7) = {X0_c + R_c, X1_c, 0., SizeInter};
Point(8) = {X0_c, X1_c + R_c, 0., SizeInter};
Point(9) = {X0_c - R_c, X1_c, 0., SizeInter};

// ... Outer Boundary
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Curve Loop(1) = {1, 2, 3, 4};

// ...  Circle segments for the inclusion
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
Curve Loop(2) = {5, 6, 7, 8};

// ... Inclusion circle:
Plane Surface(1) = {2};
// ... Rest of the domain (with Curve Loop(2) as a hole):
Plane Surface(2) = {1, 2};

// Label the outer boundaries:
Physical Curve("Top") = {3};
Physical Curve("Left") = {4};
Physical Curve("Right") = {2};
Physical Curve("Bottom") = {1};

// Label the subdomain:
Physical Surface("Inclusion") = {1};
Physical Surface("Rest") = {2};




