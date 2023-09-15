// Cartesian shock tube with two chambers

// Parameters
lx = 1;
lxp = 0.5;
ly = 0.5;
theta = 45;
NxL = 51;
NxR = 51;
Ny = 21;
dx = 0.01; // Useful for triangles only

// Mesh
// x' = x Cos(theta) - y Sin(theta)
// y' = x Sin(theta) + y Cos(theta)
theta = theta * Pi / 180;
Point(1) = {0, 0, 0, dx};
// Point(2) = {lxp, 0, 0};
Point(2) = {lxp * Cos(theta), lxp * Sin(theta), 0, dx};
// Point(3) = {lx, 0, 0};
Point(3) = {lx * Cos(theta), lx * Sin(theta), 0, dx};
// Point(4) = {lx, ly, 0};
Point(4) = {lx * Cos(theta) - ly * Sin(theta), lx * Sin(theta) + ly * Cos(theta), 0, dx};
// Point(5) = {lxp, ly, 0};
Point(5) = {lxp * Cos(theta) - ly * Sin(theta), lxp * Sin(theta) + ly * Cos(theta), 0, dx};
// Point(6) = {0, ly, 0};
Point(6) = {- ly * Sin(theta), ly * Cos(theta), 0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line(7) = {2, 5};

// Comment if you want unstructured triangles
Transfinite Line{1, 5}  = NxL;
Transfinite Line{2, 4} = NxR;
Transfinite Line{3, 6, 7}  = Ny;

Line Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};

Line Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

// Uncomment if you want unstructured triangles
Transfinite Surface{1};
Recombine Surface{1};
Transfinite Surface{2};
Recombine Surface{2};

Physical Line(1) = {3, 6};
Physical Line(2) = {1, 2, 4, 5};

Physical Surface(10) = {1}; // Left chamber
Physical Surface(11) = {2}; // Right chamber