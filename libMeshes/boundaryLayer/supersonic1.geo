// Gmsh project created on Wed Apr 06 17:26:11 2016
nbx = 31; nby = 31;
Lx = 1.; Ly = 0.6;

//Def geometrie
Point(1) = {0, 0, 0};
Point(2) = {Lx, 0, 0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0, Ly, 0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(8) = {1, 2, 3, 4};
Plane Surface(9) = {8};

//Maillage
Transfinite Line {1,-3} = nbx Using Progression 1.03;
Transfinite Line {2,-4} = nby Using Progression 1.08;
Transfinite Surface  {9};
Recombine Surface {9};

//Physical groups
Physical Line(1) = {3,4}; //Inj
Physical Line(2) = {2}; //Sortie
Physical Line(3) = {1}; //mur
surfaceFluide = 10;
Physical Surface(surfaceFluide) = {9};