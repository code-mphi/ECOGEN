cl = 0.05;
Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Line(1) = {4, 3};
Line(2) = {2, 3};
Line(3) = {1, 2};
Line(4) = {4, 1};

Physical Line(1) = {1, 3, 2, 4};
Line Loop(6) = {4, 3, 2, -1};
Plane Surface(7) = {6};
Physical Surface(8) = {7};
