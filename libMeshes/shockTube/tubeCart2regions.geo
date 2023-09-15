// Cartesian shock tube with two chambers

// Parameters
lx = 1;
lxp = 0.5;
ly = 0.5;
NxL = 51;
NxR = 51;
Ny = 51;

// Mesh
Point(1) = {0, 0, 0};
Point(2) = {lxp, 0, 0};
Point(3) = {lx, 0, 0};
Point(4) = {lx, ly, 0};
Point(5) = {lxp, ly, 0};
Point(6) = {0, ly, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line(7) = {2, 5};

Transfinite Line{1, 5}  = NxL;
Transfinite Line{2, 4} = NxR;
Transfinite Line{3, 6, 7}  = Ny;

Line Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};

Line Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Transfinite Surface{1};
Recombine Surface{1};

Transfinite Surface{2};
Recombine Surface{2};

Physical Line(1) = {3, 6};
Physical Line(2) = {1, 2, 4, 5};

Physical Surface(10) = {1}; // Left chamber
Physical Surface(11) = {2}; // Right chamber