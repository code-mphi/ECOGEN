// --- 2D structured disc with O-grid and two fluid domains ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Mesh parameters ---
NrExt = 15;  // Nb of cell along the radius on external zone
NrInt = 25;  // Nb of cell along the radius on internal zone
Ntheta = 10;  // Nb of cells along the azimuth

// --- Design parameters ---
rExt = 0.2;      // External radius of the disc
rInt = rExt/1.5;  // Internal radius for fluid domain interface
rMin = rInt/10.; // Minimum radius for O-grid construction

// --- Boundary conditions ---
boundCondWalls = 1;
fluidInt = 11;
fluidExt = 10;

// --- Geometry ---

// Center of all circles
Point(1) = {refX, refY, refZ};

// External circle
Point(2) = {refX+rExt, refY, refZ};
Point(3) = {refX, refY+rExt, refZ};
Point(4) = {refX-rExt, refY, refZ};
Point(5) = {refX, refY-rExt, refZ};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// Internal circle
Point(6) = {refX+rInt, refY, refZ};
Point(7) = {refX, refY+rInt, refZ};
Point(8) = {refX-rInt, refY, refZ};
Point(9) = {refX, refY-rInt, refZ};

Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

// Linking external/internal circle's
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Structured mesh definition outside rInt
Transfinite Line{1:4} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{5:8} = Ntheta; // Nb of cells for outside rInt
Transfinite Line{9:12} = NrExt; // Nb of cells along the radius

// Surface outside rInt
Line Loop(1) = {5, -10, -1, 9}; Plane Surface(1) = {1}; 
Line Loop(2) = {6, -11, -2, 10}; Plane Surface(2) = {2};
Line Loop(3) = {7, -12, -3, 11}; Plane Surface(3) = {3};
Line Loop(4) = {8, -9, -4, 12}; Plane Surface(4) = {4};

// Structured surface outside rInt
Transfinite Surface{1:4};
Recombine Surface{1:4};

// Structured mesh definition inside rInt
//O-grid
Point(10) = {refX+rMin, refY, refZ};
Point(11) = {refX, refY+rMin, refZ};
Point(12) = {refX-rMin, refY, refZ};
Point(13) = {refX, refY-rMin, refZ};
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {12, 13};
Line(16) = {13, 10};
// Linking internal circle and O-grid
Line(17) = {6, 10};
Line(18) = {7, 11};
Line(19) = {8, 12};
Line(20) = {9, 13};
// Structured mesh
Transfinite Line{5:8} = Ntheta;   // Nb of cells for inside rInt
Transfinite Line{13:16} = Ntheta; // Nb of cells for O-grid
Transfinite Line{17:20} = NrInt;  // Nb of cells along the radius
// Surface definition between rInt and rMin
Line Loop(5) = {5, 18, -13, -17}; Plane Surface(5) = {5};
Line Loop(6) = {6, 19, -14, -18}; Plane Surface(6) = {6};
Line Loop(7) = {7, 20, -15, -19}; Plane Surface(7) = {7};
Line Loop(8) = {8, 17, -16, -20}; Plane Surface(8) = {8};
// Structured surface between rInt and rMin
Transfinite Surface{5:8};
Recombine Surface{5:8};
// Structured mesh inside rMin
Line Loop(9) = {13, 14, 15, 16}; Plane Surface(9) = {9};
Transfinite Surface{9};
Recombine Surface{9};

// --- Physical boundaries/volume ---
Physical Line(boundCondWalls) = {1, 2, 3, 4};
Physical Surface(fluidInt) = {5, 6, 7, 8, 9};
Physical Surface(fluidExt) = {1, 2, 3, 4};