// Cartesian shock tube

// Parameters
lx = 1;
ly = 0.2;
Nx = 11;
Ny = 1;

// Mesh
Point(1) = {0, 0, 0};
Point(2) = {lx, 0, 0};
Point(3) = {lx, ly, 0};
Point(4) = {0, ly, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Transfinite Line{1,3}  = Nx;
Transfinite Line{2,4}  = Ny;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface{1};
Recombine Surface{1};

Physical Line(1) = {1,3};
Physical Line(2) = {4};
Physical Line(3) = {2};

Physical Surface(10) = {1};