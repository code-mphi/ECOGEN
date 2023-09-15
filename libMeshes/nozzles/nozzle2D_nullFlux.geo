// 2D mesh of a half nozzle

//    rE _                                          
//    rS _ |\__  _  _  _  _  _   
//         |   \__          ___/|
//         |      \__   ___/    |
//    rC _ | _  _  _ \ /        |
//     0 _ |__________|_________|
//
//         <---------><-------->
//             lConv     lDiv    

// Geometry details

// rE _ 4                                         
//      |\__
//      |   \__                                    
// rS _ | _  _ \__ (3)_  _  _   _  _  _  _   6
//      |         \__                (7)____/|
//      |(4)         \__           ____/     |
//      |               \__   ____/          |
// rC _ | _  _  _  _  _  _ \3/               |
//      |                                    |(6)
//      |       "1"      (2)|       "2"      |
//      |                                    |
//  0 _ |___________________|________________|
//      1        (1)        2      (5)       5 
//      <-------------------><------------->
//               lConv             lDiv 

//----------------- DATA -----------------

// Geo parameters
rE = 0.5;
rC = 0.4;
rS = 0.5;
lConv = 0.5;
lDiv = 0.5;

// Number of cells x-dir.
nbxLConv = 20;
nbxLDiv  = 20;

// Number of cells y-dir.
nby = 1;

// Boundary condition numbers
bcAxis = 1;
bcWall = 2;
bcInlet = 3;
bcOutlet = 4;
surfFluid = 10;

//----------------- GEO -----------------

// Geometry convergent
Point(1) = {0, 0, 0};
Point(2) = {lConv, 0, 0};
Point(3) = {lConv, rC, 0};
Point(4) = {0, rE, 0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Geometry divergent
Point(5) = {lConv+lDiv, 0, 0};
Point(6) = {lConv+lDiv, rS, 0};
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,3};
Line Loop(2) = {5,6,7,-2};
Plane Surface(2) = {2};

// Structured mesh
Transfinite Line {1,3} = nbxLConv;
Transfinite Line {5,7} = nbxLDiv;
Transfinite Line {2,4,6} = nby;
Transfinite Surface "*";
Recombine Surface "*";

// Boundary conditions and fluid domain
Physical Line(bcAxis) = {1,5};
Physical Line(bcWall) = {3,7};
Physical Line(bcInlet) = {4};
Physical Line(bcOutlet) = {6};
Physical Surface(surfFluid) = {1,2};
