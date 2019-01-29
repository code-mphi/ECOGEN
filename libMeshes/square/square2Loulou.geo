// Gmsh project created on Tue Nov 18 13:54:37 2014
dx = 0.01;

Point(1) = {0, 0, 0, dx};
Point(2) = {2, 0, 0, dx};
Point(3) = {2, 2, 0, dx};
Point(4) = {0, 2, 0, dx};
Line(1) = {4, 3};
Line(2) = {2, 3};
Line(3) = {2, 1};
Line(4) = {1, 4};
Point(5) = {1, 1, 0, dx};
Point(6) = {1.2, 1, 0, dx};
Point(7) = {1.2, 1.2, 0, dx};
Point(8) = {1, 1.2, 0, dx};
Point(9) = {0.8, 1, 0, dx};
Point(10) = {1, 0.8, 0, dx};
Delete {
  Point{7};
}
Circle(5) = {8, 5, 6};
Circle(6) = {6, 5, 10};
Circle(7) = {10, 5, 9};
Circle(8) = {9, 5, 8};
Line Loop(9) = {1, -2, 3, 4};
Line Loop(10) = {8, 5, 6, 7};
Plane Surface(11) = {9, 10};
Physical Surface(10) = {11};
Physical Line(1) = {1, 4, 3, 2};
Physical Line(2) = {6, 5, 8, 7};
