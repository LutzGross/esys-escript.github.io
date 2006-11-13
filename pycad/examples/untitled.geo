Point(1) = {0,0,0,0.1};
Point(2) = {1.6,0.8,0,0.1};
Line(1) = {1,2};
Extrude {0,0,1.} {
  Line{1};
}
Extrude {1,0,1.} {
  Surface{5};
}
Extrude {1,0,1.} {
  Surface{22,18,5,14,27};
}
