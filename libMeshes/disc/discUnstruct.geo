// --- 2D unstructured circle ---

// --- Reference frame ---
refX = 0;
refY = 0;
refZ = 0;

// --- Geometry ---
r = 0.2;

// --- Boundary conditions ---
bcWall = 1;
fluid = 10;

// --- Mesh parameter ---
dx = 0.01;

// --- Geometry ---

// Circle
Point(1) = {0, 0, 0, dx};
Point(2) = {-r, 0, 0, dx};
Point(3) = {0, r, 0, dx};
Point(4) = {r, 0, 0, dx};
Point(5) = {0, -r, 0, dx};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(1) = {1,2,3,4};

// Physical boundary/volume
Plane Surface(1) = {1}; 
Physical Line(bcWall) = {1,2,3,4};
Physical Surface(fluid) = {1};