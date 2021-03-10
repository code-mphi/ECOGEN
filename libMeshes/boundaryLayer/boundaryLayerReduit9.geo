// Gmsh project created on Wed Apr 06 17:26:11 2016
nbx1 = 16; nbx2 = 46; nby1 = 46; nby2 = 31;
Lx1 = 0.05; Lx2 = 0.20; Ly1 = 0.002; Ly2 = 0.098;

//Def geometrie
Point(1) = {0, 0, 0};
Point(2) = {Lx1, 0, 0};
Point(3) = {Lx1+Lx2, 0, 0};
Point(4) = {Lx1+Lx2, Ly1, 0};
Point(5) = {Lx1+Lx2, Ly1+Ly2, 0};
Point(6) = {Lx1, Ly1+Ly2, 0};
Point(7) = {0, Ly1+Ly2, 0};
Point(8) = {0, Ly1, 0};
Point(9) = {Lx1, Ly1, 0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};
Line(9) = {2,9};
Line(10) = {9,6};
Line(11) = {8,9};
Line(12) = {9,4};

Line Loop(13) = {1,9,-11,8};
Plane Surface(1) = {13};
Line Loop(14) = {2,3,-12,-9};
Plane Surface(2) = {14};
Line Loop(15) = {12,4,5,-10};
Plane Surface(3) = {15};
Line Loop(16) = {11,10,6,7};
Plane Surface(4) = {16};

//Maillage
Transfinite Line {1,11,-6} = nbx1;
Transfinite Line {2,12,-5} = nbx2;
Transfinite Line {-8,9,3} = nby1;
Transfinite Line {-7,10,4} = nby2 Using Progression 1.17;
Transfinite Surface  {1,2,3,4};
Recombine Surface {1,2,3,4};

//Physical groups
Physical Line(1) = {7,8}; //Inj
Physical Line(2) = {3,4,5,6}; //Sortie
Physical Line(3) = {2}; //mur
Physical Line(4) = {1}; //Abs
surfaceFluide = 10;
Physical Surface(surfaceFluide) = {1,2,3,4};
