// --- 3D structured cylinder with O-grid ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Mesh parameters ---
Nr = 10;     // Nb of cell along the radius
Ntheta = 5;  // Nb of cells along the azimuth
Nz = 25;     // Nb of cells along cylinder height

// Design parameters
r = 0.2;
rOgrid = r/3.;
h = 1.;

// --- Boundary conditions ---
boundCondInflow = 1;
boundCondOutflow = 2;
boundCondWalls = 3;
fluidVolume = 10;

// --- Geometry ---

// Inflow circle center
posZInflow = refZ;
Point(1) = {refX, refY, posZInflow};

// Inflow circle
Point(2) = {refX+r, refY, posZInflow};
Point(3) = {refX, refY+r, posZInflow};
Point(4) = {refX-r, refY, posZInflow};
Point(5) = {refX, refY-r, posZInflow};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// O-grid 
Point(6) = {refX+rOgrid, refY, posZInflow};
Point(7) = {refX, refY+rOgrid, posZInflow};
Point(8) = {refX-rOgrid, refY, posZInflow};
Point(9) = {refX, refY-rOgrid, posZInflow};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

// Linking inflow circle and O-grid
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Structured mesh definition outside O-grid
Transfinite Line{1:4} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{5:8} = Ntheta; // Nb of cells for each O-grid side
Transfinite Line{9:12} = Nr;    // Nb of cells along the radius

// Cross surface inflow outside O-grid
Line Loop(1) = {5, -10, -1, 9}; Plane Surface(1) = {1}; 
Line Loop(2) = {6, -11, -2, 10}; Plane Surface(2) = {2};
Line Loop(3) = {7, -12, -3, 11}; Plane Surface(3) = {3};
Line Loop(4) = {8, -9, -4, 12}; Plane Surface(4) = {4};

// Structured inflow surface outside O-grid
Transfinite Surface{1:4};
Recombine Surface{1:4};

// Structured mesh definition inside O-grid
Line Loop(5) = {5, 6, 7, 8}; Plane Surface(5) = {5};
Transfinite Surface{5};
Recombine Surface{5};

// Extend in z-direction to get cylinder
Extrude{0,0,h} { Surface{1:5}; Layers{Nz}; Recombine; }

// Physical boundaries/volume
Physical Surface(boundCondInflow) = {1,2,3,4,5}; 
Physical Surface(boundCondOutflow) = {34,56,78,100,122};
Physical Surface(boundCondWalls) = {29,51,73,95};

Surface Loop(6) = {1,2,3,4,5, 29,51,95,73, 34,56,78,100,122};
Volume(10) = {6};
Physical Volume(fluidVolume) = {10};