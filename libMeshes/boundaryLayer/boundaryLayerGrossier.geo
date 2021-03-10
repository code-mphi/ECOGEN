// Gmsh project created on Wed Apr 06 17:26:11 2016
nbx1 = 15; nbx2 = 35; nby = 50;
Lx1 = 0.15; Lx2 = 0.35; Ly = 0.3;

dx1 = Lx1/nbx1; dx2 = Lx2/nbx2;
dy = Ly/nby;

//Def geometrie
Point(1) = {0, 0, 0};
Point(2) = {0.5, 0, 0};
Point(3) = {0.5, 0.3, 0};
Point(4) = {0, 0.3, 0};
Point(5) = {0.15, 0, 0};
Point(6) = {0.15, 0.3, 0};
Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,6};
Line(5) = {6,4};
Line(6) = {4,1};
Line(7) = {5,6};

Line Loop(8) = {6, 1, 7, 5};
Plane Surface(9) = {8};
Line Loop(10) = {7, -4, -3, -2};
Plane Surface(11) = {10};

//Maillage
Transfinite Line {1,5} = nbx1;
Transfinite Line {2,4} = nbx2;
Transfinite Line {-6,7,3} = nby Using Progression 1.04;
Transfinite Surface  {9,11};
Recombine Surface {9,11};

//Physical groups
Physical Line(1) = {6}; //Inj
Physical Line(2) = {3, 4, 5}; //Sortie (ou abs)
Physical Line(3) = {2}; //mur
Physical Line(4) = {1}; //Abs
surfaceFluide = 10;
Physical Surface(surfaceFluide) = {9, 11};