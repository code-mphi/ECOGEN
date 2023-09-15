// 2D mesh of a smooth cross section variation nozzle

// rE _ 1                                         
//      |\__
//      |   \__                                    
// rS _ | _  _ \__ (4)_  _  _   _  _  _  _   3
//      |         \__                (5)____/|
//      |            \__           ____/     |
//      |               \__   ____/          |
// rC _ | _  _  _  _  _  _ \2/               |
//      |                                    |
//      |       "1"      (2)|       "2"      |
//      |                                    |
//  0 _ |<-------lConv----->|<-----lDiv----->|
//      |                                    |
//      |                   |                |
//      | (1)                 (2)            | (3)
//      |                   |                |
//      |                __/5\____           |
//      |             __/         \____      |
//      |          __/                 \____ |
//      |       __/                  (7)    \|
//      |    __/   (6)                       6   
//      | __/                                   
//      |/                                    
//      4

//----------------- DATA -----------------

// Geo parameters
rE = 0.14657/2.; 
rC = 0.06406/2.; 
rS = 0.14657/2.;
lConv = 0.5; 
lDiv = 0.5;

// Cells x-direction
nbxLConv = 51;
nbxLDiv  = 51;

// Cells y-direction
nby = 1;

// Boundary condition numbers
bcInlet = 1;
bcOutlet = 2;
bcWall = 3;
surfFluid = 10;

//----------------- GEO -----------------

// Upper body
Point(1) = {0, rE, 0};
Point(2) = {lConv, rC, 0};
Point(3) = {lConv+lDiv, rS, 0};

// Lower body
Point(4) = {0, -rE, 0};
Point(5) = {lConv, -rC, 0};
Point(6) = {lConv+lDiv, -rS, 0};

// Vertical links
Line(1) = {1,4};
Line(2) = {2,5};
Line(3) = {3,6};

// Close upper body
Line(4) = {1,2};
Line(5) = {2,3};

// Close lower body
Line(6) = {4,5};
Line(7) = {5,6};

// Surface convergent
Line Loop(1) = {1, 6, -2, -4};
Plane Surface (1) = {1};

// Surface divergent
Line Loop(2) = {2, 7, -3, -5};
Plane Surface (2) = {2};

// Structured mesh
Transfinite Line {4,6} = nbxLConv;
Transfinite Line {5,7} = nbxLDiv;
Transfinite Line {1,2,3} = nby;
Transfinite Surface "*";
Recombine Surface "*";

// Boundary conditions and fluid domain
Physical Line(bcInlet) = {1};
Physical Line(bcOutlet) = {3};
Physical Line(bcWall) = {4,5,6,7};
Physical Surface(surfFluid) = {1,2};