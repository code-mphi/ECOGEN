// 2D basic fan

// --- Parameters ---
// Geometry
r = 0.5;
rAxe = 0.1;
rExt = 0.45;
eBlade = 0.05;
nbBlades = 5;
// Mesh
dx = 0.025;
// Boundary conditions
boundCondExt = 1;
boundCondBlades = 2;
fluidSurface = 10;

// --- Geometry ---
// Angles
theta = 2. * Pi / nbBlades; // Between 2 blades center
phi1 = Asin(0.5*eBlade/rAxe); // Borders on center
phi2 = Asin(0.5*eBlade/rExt); // Borders at the tip of the blade

// Boundaries
Point(1) = {0, 0, 0, dx};
Point(2) = {-r, 0, 0, dx};
Point(3) = {0, r, 0, dx};
Point(4) = {r, 0, 0, dx};
Point(5) = {0, -r, 0, dx};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

//Blade after blade
angle = 0;
For k In {0:nbBlades-1}
	Point(6+k*4) = {rAxe*Sin(angle-phi1), rAxe*Cos(angle-phi1), 0, dx};
	Point(7+k*4) = {rAxe*Sin(angle+phi1), rAxe*Cos(angle+phi1), 0, dx};
	Point(8+k*4) = {rExt*Sin(angle-phi2), rExt*Cos(angle-phi2), 0, dx};
	Point(9+k*4) = {rExt*Sin(angle+phi2), rExt*Cos(angle+phi2), 0, dx};
	angle = angle + theta;
EndFor

For k In {0:nbBlades-2}
	Line(5+k*4) = {6+k*4, 8+k*4};
	Line(6+k*4) = {8+k*4, 9+k*4};
 	Line(7+k*4) = {9+k*4, 7+k*4};
 	Circle(8+k*4) = {7+k*4, 1, 6+(k+1)*4};
EndFor

k = nbBlades-1;
Line(5+k*4) = {6+k*4, 8+k*4};
Line(6+k*4) = {8+k*4, 9+k*4};
Line(7+k*4) = {9+k*4, 7+k*4};
Circle(8+k*4) = {7+k*4, 1, 6};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5:8+(nbBlades-1)*4};
Plane Surface(1) = {1,2};

Physical Surface(fluidSurface) = {1};
Physical Line(boundCondExt) = {1,2,3,4};
Physical Line(boundCondBlades) = {5:8+(nbBlades-1)*4};
