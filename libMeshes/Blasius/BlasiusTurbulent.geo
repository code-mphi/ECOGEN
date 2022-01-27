// Non-slipping wall for Blasius LEIS

// --- Geometry parameters ---
lx = 0.5;
ly = 0.3;
xw = 0.02;

// --- Mesh parameters ---
Nxl = 9;
Nxr = 193;
Ny = 91;
qy = 1.045; // Geometric progression y-direction

// --- Boundaries ---
bc_left = 1;
bc_right = 2;
bc_top = 3;
bc_down_left = 4;
bc_down_right = 5;
fluid = 10;

// --- Geometry ---
Point(1) = {0., 0., 0.};
Point(2) = {xw, 0., 0.};
Point(3) = {lx, 0., 0.};
Point(4) = {lx, ly, 0.};
Point(5) = {xw, ly, 0.};
Point(6) = {0., ly, 0.};

For i In {1:5}
    Line(i) = {i, i+1};
EndFor
Line(6) = {6, 1}; // Close outter rectangle
Line(7) = {2, 5}; // Close separation left/right rectangle

Transfinite Line{1, 5} = Nxl;
Transfinite Line{2, 4} = Nxr;
Transfinite Line{3, -6, 7} = Ny Using Progression qy;

Line Loop(1) = {2, 3, 4, -7};
Line Loop(2) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Surface{1, 2}; // No recombine to keep triangles
Recombine Surface{1, 2};

Physical Line(bc_left) = {6};
Physical Line(bc_right) = {3};
Physical Line(bc_top) = {4, 5};
Physical Line(bc_down_left) = {1};
Physical Line(bc_down_right) = {2};
Physical Surface(fluid) = {1, 2};