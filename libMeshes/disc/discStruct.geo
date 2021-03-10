// --- 2D structured circle with O-grid ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Mesh parameters ---
Nr = 10;     // Nb of cell along the radius
Ntheta = 5;  // Nb of cells along the azimuth

// --- Geometry parameters ---
r = 0.2;
rOgrid = r/3.;

// --- Boundary conditions ---
boundCondWalls = 1;
fluid = 10;

// --- Geometry ---

// Circles center
Point(1) = {refX, refY, refZ};

// External circle
Point(2) = {refX+r, refY, refZ};
Point(3) = {refX, refY+r, refZ};
Point(4) = {refX-r, refY, refZ};
Point(5) = {refX, refY-r, refZ};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// O-grid 
Point(6) = {refX+rOgrid, refY, refZ};
Point(7) = {refX, refY+rOgrid, refZ};
Point(8) = {refX-rOgrid, refY, refZ};
Point(9) = {refX, refY-rOgrid, refZ};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

// Linking external circle and O-grid
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Structured mesh definition outside O-grid
Transfinite Line{1:4} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{5:8} = Ntheta; // Nb of cells for each O-grid side
Transfinite Line{9:12} = Nr;    // Nb of cells along the radius

// Surface outside O-grid
Line Loop(1) = {5, -10, -1, 9}; Plane Surface(1) = {1}; 
Line Loop(2) = {6, -11, -2, 10}; Plane Surface(2) = {2};
Line Loop(3) = {7, -12, -3, 11}; Plane Surface(3) = {3};
Line Loop(4) = {8, -9, -4, 12}; Plane Surface(4) = {4};

// Structured surface outside O-grid
Transfinite Surface{1:4};
Recombine Surface{1:4};

// Structured mesh definition inside O-grid
Line Loop(5) = {5, 6, 7, 8}; Plane Surface(5) = {5};
Transfinite Surface{5};
Recombine Surface{5};

// Physical boundary/surface
Physical Line(boundCondWalls) = {1, 2, 3, 4};
Physical Surface(fluid) = {1, 2, 3, 4, 5};