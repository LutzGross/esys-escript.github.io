lc_border = 0.4;
lc_surf = 0.4;

N=8;
DTheta = 2*Pi/N;
base_z = 0.0;
top_z = 4.0;
surf_top_z = 7.0;
Theta = 0;
rad_conduit = 1.5;
rad_surf = 3.0;

CenterPoint1 = newp;
Point(CenterPoint1) = { 0.0,0.0,base_z,lc_border};

// Conduit bottom Circle
For I In {1:N}
   PointList[I] = newp;
   Theta+=DTheta;
   Point( PointList[ I] ) = { rad_conduit*Cos(Theta)  
,rad_conduit*Sin(Theta), base_z, lc_border};
EndFor

// Conduit Top Circle
CenterPoint2 = newp;
Point(CenterPoint2) = {0.0,0.0,top_z,lc_border};
Theta = 0.0;
For I In {1:N}
   PointList2[I] = newp;
   Theta+=DTheta;
   Point( PointList2[ I] ) = { rad_conduit*Cos(Theta)  
,rad_conduit*Sin(Theta), top_z, lc_border};
EndFor

// Bottom Surface Lines
For I In {1:N}
   LineList[I] = newl;
   Line( LineList[I] ) = { PointList[I], CenterPoint1 };
EndFor
For I In {1:N}
   CircleList[I] = newc;
   Circle( CircleList[I] ) = { PointList[I], CenterPoint1,  
PointList[(I%N) + 1] };
EndFor


// Top Surface Lines
//For I In {1:N}
//  LineList2[I] = newl;
//  Line( LineList2[I] ) = { PointList2[I], CenterPoint2 };
//EndFor
For I In {1:N}
   CircleList2[I] = newc;
   Circle( CircleList2[I] ) = { PointList2[I], CenterPoint2,  
PointList2[(I%N) + 1] };
EndFor



// Lines Connecting Top and Bottom Surfaces
For I In {1:N}
   Line(newl) = { PointList[I], PointList2[I] };
EndFor



// ******   DOME ****************
// ******************************

// Surf Bot
Theta = 0.0;
For I In {1:N}
   SurfPoint1[I] =  newp;
   Theta += DTheta;
   Point( SurfPoint1[I] ) = { rad_surf*Cos(Theta)  
,rad_surf*Sin(Theta), top_z, lc_surf };
EndFor



// Surf Top
CenterPoint3 = newp;
Point(CenterPoint3) = {0.0,0.0,surf_top_z,lc_border};
Theta = 0.0;
For I In {1:N}
   SurfPoint2[I] = newp;
   Theta+=DTheta;
   Point( SurfPoint2[ I] ) = { rad_surf*Cos(Theta)  
,rad_surf*Sin(Theta), surf_top_z, lc_surf};
EndFor

// Lines Connecting Top and bottom
For I In {1:N}
   Line(newl) = { SurfPoint1[I], SurfPoint2[I] };
EndFor

// Lines connecting outer to inner circle on surface

For I In {1:N}
   Line(newl) = { PointList2[I], SurfPoint1[I]};
EndFor



// Bottom Surface Lines
For I In {1:N}
   CircleList4[I] = newc;
   Circle(CircleList4[I] ) = { SurfPoint1[I], CenterPoint2,  
SurfPoint1[(I%N) + 1] };
EndFor



// Top Surface Lid
For I In {1:N}
   LineList3[I] = newl;
   Line( LineList3[I] ) = { SurfPoint2[I], CenterPoint3 };
EndFor
For I In {1:N}
   CircleList3[I] = newc;
   Circle( CircleList3[I] ) = { SurfPoint2[I], CenterPoint3,  
SurfPoint2[(I%N) + 1] };
EndFor

// ** BEGIN **
// conduit Bottom Surface
Line Loop(73) = {4,-5,-12};
Plane Surface(74) = {73};
Line Loop(75) = {3,-4,-11};
Plane Surface(76) = {75};
Line Loop(77) = {10,3,-2};
Plane Surface(78) = {77};
Line Loop(79) = {9,2,-1};
Plane Surface(80) = {79};
Line Loop(81) = {16,1,-8};
Plane Surface(82) = {81};
Line Loop(83) = {15,8,-7};
Plane Surface(84) = {83};
Line Loop(85) = {14,7,-6};
Plane Surface(86) = {85};
Line Loop(87) = {5,-6,-13};
Plane Surface(88) = {87};


Physical Surface(153) = {82,80,78,76,74,88,86,84};
//  ** END **

// ** BEGIN **
// Volcano Top Surface
Line Loop(89) = {54,-47,-22,46};
Plane Surface(90) = {89};
Line Loop(91) = {48,-55,-47,23};
Plane Surface(92) = {91};
Line Loop(93) = {24,41,-56,-48};
Plane Surface(94) = {93};
Line Loop(95) = {49,-42,-17,41};
Plane Surface(96) = {95};
Line Loop(97) = {50,-43,-18,42};
Plane Surface(98) = {97};
Line Loop(99) = {51,-44,-19,43};
Plane Surface(100) = {99};
Line Loop(101) = {45,-52,-44,20};
Plane Surface(102) = {101};
Line Loop(103) = {53,-46,-21,45};
Plane Surface(104) = {103};


Physical Surface(154) = {90,92,94,96,104,102,100,116,98};
// ** END **

// ** BEGIN **
// Volcano Free Surface Top
Line Loop(105) = {70,63,-62};
Plane Surface(106) = {105};
Line Loop(107) = {62,-61,69};
Plane Surface(108) = {107};
Line Loop(109) = {61,-60,68};
Plane Surface(110) = {109};
Line Loop(111) = {60,-59,67};
Plane Surface(112) = {111};
Line Loop(113) = {58,-59,-66};
Plane Surface(114) = {113};
Line Loop(115) = {57,-58,-65};
Plane Surface(116) = {115};
Line Loop(117) = {72,57,-64};
Plane Surface(118) = {117};
Line Loop(119) = {63,-64,-71};
Plane Surface(120) = {119};


Physical Surface(155) = {106,110,112,114,116,118,120,108};
// ** END **


// ** BEGIN **
// Conduit Outer Surface
Line Loop(121) = {13,30,-21,-29};
Ruled Surface(122) = {121};
Line Loop(123) = {29,-20,-28,12};
Ruled Surface(124) = {123};
Line Loop(125) = {28,-19,-27,11};
Ruled Surface(126) = {125};
Line Loop(127) = {27,-18,-26,10};
Ruled Surface(128) = {127};
Line Loop(129) = {26,-17,-25,9};
Ruled Surface(130) = {129};
Line Loop(131) = {25,-24,-32,16};
Ruled Surface(132) = {131};
Line Loop(133) = {32,-23,-31,15};
Ruled Surface(134) = {133};
Line Loop(135) = {30,22,-31,-14};
Ruled Surface(136) = {135};


Physical Surface(156) = {134,136,122,126,128,130,132,124};
// ** END **


// ** BEGIN **
// Outer Free Surface Cylinder
Line Loop(137) = {53,38,-69,-37};
Ruled Surface(138) = {137};
Line Loop(139) = {52,37,-68,-36};
Ruled Surface(140) = {139};
Line Loop(141) = {51,36,-67,-35};
Ruled Surface(142) = {141};
Line Loop(143) = {50,35,-66,-34};
Ruled Surface(144) = {143};
Line Loop(145) = {49,34,-65,-33};
Ruled Surface(146) = {145};
Line Loop(147) = {56,33,-72,-40};
Ruled Surface(148) = {147};
Line Loop(149) = {55,40,-71,-39};
Ruled Surface(150) = {149};
Line Loop(151) = {39,-70,-38,54};
Ruled Surface(152) = {151};

Physical Surface(157) = {138,148,144,142,140,150,146,152};
// ** End **


Surface Loop(158) =  
{134,132,130,128,126,124,122,88,74,76,-78,-80,-82,-84,-86,-136,-90,152,150,92,94,-96,146,144,-98,-100,142,140,102,-104,138,108,106,-120,118,148,-116,-114,112,110};

Volume(159) = {158};
Physical Volume(160)  = { 159};
