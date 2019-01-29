// Gmsh project created on Fri Sep 12 15:39:15 2014
k=0.05;
Point(1) = {0, 0, 0, k};
Point(2) = {2, 0, 0, k};
Point(3) = {3, 2, 0, k};
Point(4) = {1, 2, 0, k};
Point(5) = {1.5, 1, 0, k};

Spline(2) = {4, 5, 1, 2};
Spline(3) = {2, 3, 4};
Line Loop(3) = {2,3};
Plane Surface(1) = {3};

Extrude {0, 0, 1} {
  Surface{1};
}
Physical Surface(1) = {14, 15};
Physical Surface(2) = {1, 10};
Physical Volume(17) = {1};

//Recombine Surface {14, 15, 10, 1};
