// Gmsh project created on Tue Sep 09 13:41:01 2014
dx = 0.05;
dxMarche = 0.025;

Point(1) = {0, 0, 0, dx};
Point(2) = {2, 0, 0, dx};
Point(3) = {2, 0.3, 0, dx};
Point(4) = {0, 0.3, 0, dx};
Point(5) = {1.5, 0, 0, dxMarche};
Point(6) = {1.52, 0, 0, dxMarche};
Point(7) = {1.5, 0.05, 0, dxMarche};
Point(8) = {1.52, 0.05, 0, dxMarche};


Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 6};
Line(4) = {1, 4};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
Line(8) = {1, 5};
Line Loop(9) = {4, 1, 2, 3, -7, -6, -5, -8};
Plane Surface(10) = {9};
Physical Surface(11) = {10};
Physical Line(1) = {8, 4, 1, 3, 7, 6, 5};
Physical Line(2) = {2};

Recombine Surface {10};
