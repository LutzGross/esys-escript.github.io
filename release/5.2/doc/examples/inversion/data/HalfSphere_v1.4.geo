// Half-sphere example from Rucker et al (2006), Figure 4.
// "Three-dimensional modelling and inversion of dc resistivity data incorporating topography – I. Modelling"
// Geophys. J. Int. (2006) 166, 495–505
// 
// 


// ------------------------------------------------------------------------------------------------
// DEFINITIONS
// ------------------------------------------------------------------------------------------------

// Modelling domain:
xmin_domain =-100;
xmax_domain = 100;
ymin_domain =-100; 
ymax_domain = 100;   
zmin_domain =-100; 
zmax_domain =   0;

// Sphere:
R = 2.25; // sphere radius
x_sphere = (xmin_domain + xmax_domain) / 2.0; // sphere center X
y_sphere = (ymin_domain + ymax_domain) / 2.0; // sphere center Y
z_sphere = zmax_domain; // sphere top

// Local element sizes:
cl_0 =  25.0; // domain
cl_1 =   0.5; // sphere
ep_0 =   0.1;

// Small value to add to electrode position die to overlap with prism-points:
eps  = 1.0e-2  ; 

// Electrodes:
xe_0 = -5  ; // start X-coordinate
ye_0 =  0  ; // start Y-coordinate
ze_0 =  0  ; // start Z-coordinate
numx = 21  ; // number of electrides in X-direction
numy =  1  ; // number of electrides in Y-direction
estp = 0.5 ; // step size



//-------------------------------------------------------------------------------------------------
//  DOMAIN
//-------------------------------------------------------------------------------------------------

// Points top, domain
Point(01) = { xmin_domain, ymin_domain, zmax_domain, cl_0 };
Point(02) = { xmax_domain, ymin_domain, zmax_domain, cl_0 };
Point(03) = { xmax_domain, ymax_domain, zmax_domain, cl_0 };
Point(04) = { xmin_domain, ymax_domain, zmax_domain, cl_0 };
// Points bottom, domain
Point(05) = { xmin_domain, ymin_domain, zmin_domain, cl_0 };
Point(06) = { xmax_domain, ymin_domain, zmin_domain, cl_0 };
Point(07) = { xmax_domain, ymax_domain, zmin_domain, cl_0 };
Point(08) = { xmin_domain, ymax_domain, zmin_domain, cl_0 };
 
// lines top          // lines bottom       // lines top-bottom 
Line(01)  = { 01, 02 }; Line(05)  = { 05, 06 }; Line(09)  = { 01, 05 };
Line(02)  = { 02, 03 }; Line(06)  = { 06, 07 }; Line(10)  = { 02, 06 };
Line(03)  = { 03, 04 }; Line(07)  = { 07, 08 }; Line(11)  = { 03, 07 };
Line(04)  = { 04, 01 }; Line(08)  = { 08, 05 }; Line(12)  = { 04, 08 };


//-------------------------------------------------------------------------------------------------
//  SPHERE
//-------------------------------------------------------------------------------------------------

// Points sphere
Point(10) = { x_sphere   , y_sphere   , z_sphere   , cl_1 }; // centre
Point(11) = { x_sphere-R , y_sphere   , z_sphere   , cl_1 }; // left
Point(12) = { x_sphere   , y_sphere   , z_sphere-R , cl_1 }; // bottom
Point(13) = { x_sphere+R , y_sphere   , z_sphere   , cl_1 }; // right
Point(14) = { x_sphere   , y_sphere-R , z_sphere   , cl_1 }; // front
Point(15) = { x_sphere   , y_sphere+R , z_sphere   , cl_1 }; // back

// Arcs composing the half-sphere
Circle(20) = { 11, 10, 12 }; // left-bottom
Circle(21) = { 12, 10, 13 }; // bottom-right
Circle(22) = { 11, 10, 14 }; // left-front
Circle(23) = { 14, 10, 13 }; // front-right
Circle(24) = { 13, 10, 15 }; // right-back
Circle(25) = { 15, 10, 11 }; // back-left
Circle(26) = { 14, 10, 12 }; // front-bottom
Circle(27) = { 12, 10, 15 }; // bottom-back


//-------------------------------------------------------------------------------------------------
//  SURFACES
//-------------------------------------------------------------------------------------------------

// Loop Domain
Line Loop(01) = { 01, 02, 03, 04,   // top
                 -22,-23,-24,-25 }; // hole (sphere)
Line Loop(02) = { 05, 06, 07, 08 }; // bottom
Line Loop(03) = { 09, 05,-10,-01 }; // front
Line Loop(04) = { 10, 06,-11,-02 }; // right
Line Loop(05) = { 11, 07,-12,-03 }; // back
Line Loop(06) = { 12, 08,-09,-04 }; // left 

// Surfaces Domain
Plane Surface(01) = { 01 };
Plane Surface(02) = { 02 };
Plane Surface(03) = { 03 };
Plane Surface(04) = { 04 };
Plane Surface(05) = { 05 };
Plane Surface(06) = { 06 }; 

// Loop Sphere
Line Loop(07) = { 22, 23, 24, 25}; // top 
Line Loop(08) = { 20,-26,-22    };
Line Loop(09) = { 21,-23, 26    };
Line Loop(10) = { 20, 27, 25    };
Line Loop(11) = { 21, 24,-27    };
 
// Surfaces Sphere
Ruled Surface(07) = { 07 }; // top surface
Ruled Surface(08) = { 08 };
Ruled Surface(09) = { 09 };
Ruled Surface(10) = { 10 };
Ruled Surface(11) = { 11 }; 


//-------------------------------------------------------------------------------------------------
//  VOLUMES
//-------------------------------------------------------------------------------------------------

// Domain
Surface Loop(01) = { 01,02,03,04,05,06,   // domain surface with holes:
                     08,09,10,11       }; // surface sphere (sans top, part of 'holes')
// Sphere
Surface Loop(02) = { 07,08,09,10,11 };

// Volume definitions & labels:
Volume(01) = { 01 }; Physical Volume("domain") = {01}; 
Volume(02) = { 02 }; Physical Volume("sphere") = {02}; 


//-------------------------------------------------------------------------------------------------
// POINTS ELECTRODES
//-------------------------------------------------------------------------------------------------

// Electrodes are defined as embedded points in the top surface:
x = xe_0; // X-coordinate start.
y = xe_0; // Y-coordinate start.
k = 16;   // Index value of Point. ( <--- Last point index + 1)
// Loop to setup embedded points In Surface(01):
For i In{0:numx-1}
   // X-Coordinate of point:
   x = xe_0 + i*estp;   
   For j In{0:numy-1}
      // Y-coordinate of point:
      y = ye_0 + j*estp;   
      Printf("%g %g", x, y);
      // Define the new point in the surface:
      Point(k) = {x, y, ze_0, ep_0}; Point{k} In Surface{01}; 
      // Increment point index:
      k = k + 1;
   EndFor
EndFor



//-------------------------------------------------------------------------------------------------
// FIELDS
//-------------------------------------------------------------------------------------------------

// The fields setup the level of discretisation for selected entities.
// Each electrode gets a high level of elements in the neighbourhood,
// the inner region and prism have, respectively, smaller elements and
// the domain outside has the default size as specified via 'cl_0'.

// Attractor field and Threshold around the embedded electrode points:
Field[1] = Attractor;
// List of the electrodes indices -- haven't found a way to do this via a for-loop and array!
Field[1].NodesList = {16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMax = cl_0;
Field[2].DistMin = ep_0;
Field[2].LcMin = ep_0/2;   // Set local mesh size
Field[2].DistMax = cl_0*2; // Set local Area


// Box for prism:
Field[3] = Box;
Field[3].VIn  = cl_1; // prism element size
Field[3].VOut = cl_0; // outer region element size
Field[3].XMin = x_sphere-R;
Field[3].XMax = x_sphere+R;
Field[3].YMin = y_sphere-R;
Field[3].YMax = y_sphere+R;
Field[3].ZMin = z_sphere-R;
Field[3].ZMax = z_sphere;


// Background Min field required:
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 4;

// Extend the elements sizes from the boundary inside the domain:
Mesh.CharacteristicLengthExtendFromBoundary = 1;


// End script.



