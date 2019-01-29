nbx = 50; nby = 25; nbz = 25;
Lx = 2; Ly = 1; Lz = 1;

dx = Lx/nbx;
dy = Ly/nby;
dz = Lz/nbz;

// Gmsh project created on Tue Nov 04 15:37:28 2014
Point(1) = {0, 0, 0, dz};

Extrude {0, 0, Lz} {
  Point{1}; Layers{nbz}; Recombine;
}
Extrude {0, Ly, 0} {
  Line{1}; Layers{nby}; Recombine;
}
Extrude {Lx, 0, 0} {
  Surface{5}; Layers{nbx}; Recombine;
}

Physical Surface(1) = {5};
Physical Surface(2) = {27};
Physical Surface(3) = {14};
Physical Surface(4) = {22};
Physical Surface(5) = {26};
Physical Surface(6) = {18};
Physical Volume(34) = {1};
