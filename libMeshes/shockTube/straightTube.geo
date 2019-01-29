// Gmsh project created on Tue Sep 09 13:41:01 2014
k = 0.06;
Point(1) = {0, 0, 0, k};
Point(2) = {2, 0, 0, k};
Point(3) = {2, 0.3, 0, k};
Point(4) = {0, 0.3, 0, k};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line(1) = {1, 3};
Physical Line(2) = {4, 2};
Physical Surface(9) = {6};
Recombine Surface {6};
