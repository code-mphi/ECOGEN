// 2D annular disc

// --- Geometry parameters ---
rInt = 0.1;
rExt = 0.5;

// --- Boundary conditions ---
bcWallExt = 1;
bcWallInt = 2;
fluidSurface = 10;

// --- Mesh parameter ---
dx = 0.05;

// --- Geometry ---

// Exterior circle
Point(1) = {0, 0, 0, dx};
Point(2) = {-rExt, 0, 0, dx};
Point(3) = {0, rExt, 0, dx};
Point(4) = {rExt, 0, 0, dx};
Point(5) = {0, -rExt, 0, dx};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(1) = {1,2,3,4};

// Interior circle
Point(6) = {-rInt, 0, 0, dx};
Point(7) = {0, rInt, 0, dx};
Point(8) = {rInt, 0, 0, dx};
Point(9) = {0, -rInt, 0, dx};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Line Loop(2) = {5,6,7,8};

// Physical boundaries/surface
Plane Surface(1) = {1,2}; 
Physical Line(bcWallExt) = {1,2,3,4};
Physical Line(bcWallInt) = {5,6,7,8};
Physical Surface(fluidSurface) = {1};