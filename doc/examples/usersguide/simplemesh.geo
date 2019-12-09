// dimensions and mesh size
xdim = 100.;
ydim = 200.;
zdim = 50.;
mtop = 5.;
mbase = 2.;

//Points
Point(1) = {0., 0., 0., mbase};
Point(2) = {xdim, 0., 0., mbase};
Point(3) = {0., ydim, 0., mbase};
Point(4) = {xdim, ydim, 0., mbase};
Point(5) = {0., 0., zdim, mtop};
Point(6) = {xdim, 0., zdim, mtop};
Point(7) = {0., ydim, zdim, mtop};
Point(8) = {xdim, ydim, zdim, mtop};

//Lines and surfaces
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {1, 3};
Line(4) = {2, 4};
Line(5) = {5, 6};
Line(6) = {7, 8};
Line(7) = {5, 7};
Line(8) = {6, 8};
Line(9) = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};
Line Loop(1) = {-1, 3, 2, -4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 8, -6, -7};
Plane Surface(2) = {2};
Line Loop(3) = {1, 11, -5, -9};
Plane Surface(3) = {3};
Line Loop(4) = {-2, 10, 6, -12};
Plane Surface(4) = {4};
Line Loop(5) = {-3, 9, 7, -10};
Plane Surface(5) = {5};
Line Loop(6) = {4, 12, -8, -11};
Plane Surface(6) = {6};

// domain
Surface Loop(1) = {1:6};
Volume(1) = {1};
