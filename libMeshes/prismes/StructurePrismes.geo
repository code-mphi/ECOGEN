lc = 0.1;
nbPointExtrude = 8;

//Triangle de base
Point(1)={0,0,0,lc};
Point(2)={2,0,0,lc};
Point(3)={2,1,0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,1};

//Surface de base triangle
Line Loop(1)={1,2,3};
Plane Surface(1)={1};

//Volume par extrusion par translation
Extrude {0,0,-1}{Surface{1};Layers{nbPointExtrude};Recombine;}
Physical Surface(1) = {11, 15, 19, 20, 1};
Physical Volume(22) = {1};

