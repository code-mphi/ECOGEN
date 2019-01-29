// Gmsh project created on Tue Nov 18 13:54:37 2014
dx = 0.2;

Point(1) = {0, 0, 0, dx};
Point(2) = {2, 0, 0, dx};
Point(3) = {2, 2, 0, dx};
Point(4) = {0, 2, 0, dx};
Line(1) = {4, 3};
Line(2) = {2, 3};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(5) = {1, -2, 3, 4};
Plane Surface(6) = {5};
Physical Line(1) = {1, 2, 3, 4};
Physical Surface(8) = {6};
