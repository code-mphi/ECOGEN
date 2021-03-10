// 2D annular disc

// --- Geometry parameters ---
rInt = 0.1;
rExt = 0.5;

// --- Boundary conditions ---
bcWallExt = 1;
bcWallInt = 2;
fluidSurface = 10;

// --- Mesh parameters ---
Nr = 10;     // Nb of cells along radius
Ntheta = 10; // Nb of cells along the azimuth (1/4)

// --- Geometry ---

// Exterior circle
Point(1) = {0, 0, 0};
Point(2) = {-rExt, 0, 0};
Point(3) = {0, rExt, 0};
Point(4) = {rExt, 0, 0};
Point(5) = {0, -rExt, 0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Interior circle
Point(6) = {-rInt, 0, 0};
Point(7) = {0, rInt, 0};
Point(8) = {rInt, 0, 0};
Point(9) = {0, -rInt, 0};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Linking interior/exterior circles
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Structured mesh definition
Transfinite Line{1:4} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{5:8} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{9:12} = Nr;    // Nb of cells along the radius

// Surface definition
Line Loop(1) = {1, 10, -5, -9}; Plane Surface(1) = {1}; 
Line Loop(2) = {2, 11, -6, -10}; Plane Surface(2) = {2};
Line Loop(3) = {3, 12, -7, -11}; Plane Surface(3) = {3};
Line Loop(4) = {4, 9, -8, -12}; Plane Surface(4) = {4};

// Structured surface
Transfinite Surface{1:4};
Recombine Surface{1:4};

// Physical boundaries/surface
Physical Line(bcWallExt) = {1, 2, 3, 4};
Physical Line(bcWallInt) = {5, 6, 7, 8};
Physical Surface(fluidSurface) = {1, 2, 3, 4};