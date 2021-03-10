// --- 3D unstructured cylinder ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Mesh parameter ---
dx = 0.1;

// --- Geometry parameters ---
r = 0.2;
h = 1.;

// --- Boundary conditions ---
boundCondWall = 1;
boundCondInflow = 2;
boundCondOutflow = 3;
fluidVolume = 10;

// --- Geometry ---

// Center of inflow/outflow circles
pCenterBottom = 1;
posZBottom = refZ;
Point(pCenterBottom) = {refX, refY, posZBottom, dx};

pCenterTop = 2;
posZTop = refZ + h;
Point(pCenterTop) = {refX, refY, posZTop, dx};

// Bottom circle
Point(3) = {refX + r, refY, posZBottom, dx};
Point(4) = {refX, refY + r, posZBottom, dx};
Point(5) = {refX - r, refY, posZBottom, dx};
Point(6) = {refX, refY - r, posZBottom, dx};
Circle(1) = {3, pCenterBottom, 4};
Circle(2) = {4, pCenterBottom, 5};
Circle(3) = {5, pCenterBottom, 6};
Circle(4) = {6, pCenterBottom, 3};

// Top circle
Point(7) = {refX + r, refY, posZTop, dx};
Point(8) = {refX, refY + r, posZTop, dx};
Point(9) = {refX - r, refY, posZTop, dx};
Point(10) = {refX, refY - r, posZTop, dx};
Circle(5) = {7, pCenterTop, 8};
Circle(6) = {8, pCenterTop, 9};
Circle(7) = {9, pCenterTop, 10};
Circle(8) = {10, pCenterTop, 7};

// Lateral lines linking circles
Line(9) = {3, 7};
Line(10) = {4, 8};
Line(11) = {5, 9};
Line(12) = {6, 10};

// Lateral surfaces linking circles
Line Loop(1) = {-1, 9, 5, -10};
Line Loop(2) = {-2, 10, 6, -11};
Line Loop(3) = {-3, 11, 7, -12};
Line Loop(4) = {-4, 12, 8, -9};
Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};

// Inflow/outflow surfaces
Line Loop(5) = {1, 2, 3, 4};
Line Loop(6) = {5, 6, 7, 8};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Physical boundaries/volume
Physical Surface(boundCondWall) = {1:4};
Physical Surface(boundCondInflow) = {5};
Physical Surface(boundCondOutflow) = {6};

Surface Loop(1) = {1:6};
Volume(1) = {1};
Physical Volume(fluidVolume) = {1};