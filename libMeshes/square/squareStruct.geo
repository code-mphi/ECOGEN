// Geometry parameters
lx = 1.;
ly = 1.;

// Mesh parameters
Nx = 51;
Ny = 51;
struct = 1;
structQuad = 0;

// Geometry
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

If (struct == 1)
Transfinite Line{1, 3} = Nx;
Transfinite Line{2, 4} = Ny;
EndIf

Line Loop(1) = {1:4};
Plane Surface(1) = {1};

If (struct == 1)
Transfinite Surface{1};
If (structQuad == 1)
Recombine Surface{1};
EndIf
EndIf

// Boundary conditions
Physical Line(1) = {1:4};

// Fluid
Physical Surface(10) = {1};