nbx = 100; nby = 100;
Lx = 320.e-6; Ly = 320.e-6;

dx = Lx/nbx;
dy = Ly/nby;

// Gmsh project created on Tue Nov 04 15:37:28 2014
Point(1) = {0, 0, 0, dy};

Extrude {0, Ly, 0} {
  Point{1}; Layers{nby}; Recombine;
}
Extrude {Lx, 0, 0} {
  Line{1}; Layers{nbx};
}

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(10) = {5};
