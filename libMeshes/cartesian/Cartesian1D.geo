dx = 0.001;

// Gmsh project created on Tue Nov 04 15:37:28 2014
Point(1) = {0, 0, 0, dx};
Point(2) = {2, 0, 0, dx};
Line(1) = {1, 2};
Physical Point(1) = {1};
Physical Point(2) = {2};
Physical Line(4) = {1};
