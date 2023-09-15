// --- 3D structured cylinder with O-grid ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Mesh parameters ---
Nr = 10;     // Nb of cell along the radius
Ntheta = 5;  // Nb of cells along the azimuth
Nz = 100;     // Nb of cells along cylinder height

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

// ------- 

// Outflow circle center
posZOutflow = refZ+h;
Point(10) = {refX, refY, posZOutflow};

// Outflow circle
Point(11) = {refX+r, refY, posZOutflow};
Point(12) = {refX, refY+r, posZOutflow};
Point(13) = {refX-r, refY, posZOutflow};
Point(14) = {refX, refY-r, posZOutflow};

Circle(13) = {11,10,12};
Circle(14) = {12,10,13};
Circle(15) = {13,10,14};
Circle(16) = {14,10,11};

// O-grid 
Point(15) = {refX+rOgrid, refY, posZOutflow};
Point(16) = {refX, refY+rOgrid, posZOutflow};
Point(17) = {refX-rOgrid, refY, posZOutflow};
Point(18) = {refX, refY-rOgrid, posZOutflow};

Line(17) = {15, 16};
Line(18) = {16, 17};
Line(19) = {17,18};
Line(20) = {18, 15};

// Linking outflow circle and O-grid
Line(21) = {11, 15};
Line(22) = {12, 16};
Line(23) = {13, 17};
Line(24) = {14, 18};

// Structured mesh definition outside O-grid
Transfinite Line{13:16} = Ntheta; // Nb of cells for each circle azimuth
Transfinite Line{17:20} = Ntheta; // Nb of cells for each O-grid side
Transfinite Line{21:24} = Nr;    // Nb of cells along the radius

// Cross surface outflow outside O-grid
Line Loop(6) = {17, -22, -13, 21}; Plane Surface(6) = {6}; 
Line Loop(7) = {18, -23, -14, 22}; Plane Surface(7) = {7};
Line Loop(8) = {19, -24, -15, 23}; Plane Surface(8) = {8};
Line Loop(9) = {20, -21, -16, 24}; Plane Surface(9) = {9};

// Structured outflow surface outside O-grid
Transfinite Surface{6:9};
Recombine Surface{6:9};

// Structured mesh definition inside O-grid
Line Loop(10) = {17, 18, 19, 20}; Plane Surface(10) = {10};
Transfinite Surface{10};
Recombine Surface{10};

// ------- 

// Links inflow to outflow

Line(25) = {2, 11};
Line(26) = {3, 12};
Line(27) = {4, 13};
Line(28) = {5, 14};

// Links inflow O-grid to outflow O-grid

Line(29) = {6, 15};
Line(30) = {7, 16};
Line(31) = {8, 17};
Line(32) = {9, 18};

Transfinite Line{25:32} = Nz; // Nb of cells along the cylinder length

// ------- 

// Plane surfaces outside O-grid linking inflow/outflow
Line Loop(11) = {29, -21, -25, 9}; Plane Surface(11) = {11}; 
Line Loop(12) = {30, -22, -26, 10}; Plane Surface(12) = {12};
Line Loop(13) = {31, -23, -27, 11}; Plane Surface(13) = {13};
Line Loop(14) = {32, -24, -28, 12}; Plane Surface(14) = {14};

// Structured surfaces
Transfinite Surface{11:14};
Recombine Surface{11:14};

// Plane surfaces O-grid linking inflow/outflow
Line Loop(15) = {29, 17, -30, -5}; Plane Surface(15) = {15}; 
Line Loop(16) = {30, 18, -31, -6}; Plane Surface(16) = {16};
Line Loop(17) = {31, 19, -32, -7}; Plane Surface(17) = {17};
Line Loop(18) = {32, 20, -29, -8}; Plane Surface(18) = {18};

// Structured surfaces
Transfinite Surface{15:18};
Recombine Surface{15:18};

// ------- 

// Curved surfaces cylinder linking inflow/outflow
Line Loop(19) = {25, 13, -26, -1}; Surface(19) = {19};
Line Loop(20) = {26, 14, -27, -2}; Surface(20) = {20};
Line Loop(21) = {27, 15, -28, -3}; Surface(21) = {21};
Line Loop(22) = {28, 16, -25, -4}; Surface(22) = {22};

// Structured surfaces
Transfinite Surface{19:22};
Recombine Surface{19:22};

// -------

// Volume between O-grids
Surface Loop(1) = {5, 15, 16, 17, 18, 10};
Volume(1) = {1};
Transfinite Volume{1};
Recombine Volume{1};

// Volume 1/4 cylinder
Surface Loop(2) = {1, 11, 12, 15, 19, 6};
Volume(2) = {2};
Transfinite Volume{2};
Recombine Volume{2};

// Volume 2/4 cylinder
Surface Loop(3) = {2, 12, 13, 16, 20, 7};
Volume(3) = {3};
Transfinite Volume{3};
Recombine Volume{3};

// Volume 3/4 cylinder
Surface Loop(4) = {3, 13, 14, 17, 21, 8};
Volume(4) = {4};
Transfinite Volume{4};
Recombine Volume{4};

// Volume 4/4 cylinder
Surface Loop(5) = {4, 11, 14, 18, 22, 9};
Volume(5) = {5};
Transfinite Volume{5};
Recombine Volume{5};

// Physical boundaries/volume
Physical Surface(boundCondInflow) = {1:5}; 
Physical Surface(boundCondOutflow) = {6:10};
Physical Surface(boundCondWalls) = {19:22};

Physical Volume(fluidVolume) = {1:5};