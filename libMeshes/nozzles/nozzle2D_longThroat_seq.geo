// 2D mesh of a half nozzle with long throat

//    rE _  _________________                                         
//    rS _ | _  _  _  _  _  _|\__  _  _  _  _  _  _  _  _  _   _______  
//         |                 |   \__                      ___/|       |
//         |                 |      \__               ___/    |       |
//    rC _ | _  _  _  _  _  _| _  _  _ \ ___________ /        |       |
//     0 _ |_________________|__________|___________|_________|_______|
//
//         <-----------------><---------><----------><--------><------>
//                  lE            lConv       lCol      lDiv      lS 


// Geometry details

// rE _ 4__________(3)_________3                                         
//      |                      |\__
//      |                          \__                                    
// rS _ | _  _  _  _  _  _  _  | _  _ \__ (7)_  _  _  _  _  _  _  _  _  _  _  _  _  _ 10_______(16)_______12
//      |                                \__                                      ____/                   |
//      |(4)                (2)|            \__                          (13)____/     |                  |
//      |           "1"                        \__                      ____/                             | 
// rC _ | _  _  _  _  _  _  _  | _  _  _  _  _  _ \6_______(10)_______8/               |(12)     "5"      |(15)
//      |                                                                                                 |
//      |                      |                (6)|                  |(9)             |                  |
//      |                              "2"                  "3"               "4"                         |
//  0 _ |______________________|___________________|__________________| _______________|__________________|
//      1          (1)         2         (5)       5        (8)       7       (11)     9        (14)      11
//      <----------------------><------------------><-----------------><---------------><----------------->
//                  lE                   lConv              lCol               lDiv               lS 

//----------------- DATA -----------------

// Geo parameters
rE = 0.13; rC = 0.04; rS = 0.13;
lE = 0.05; lConv = 0.2; lCol = 0.05; lDiv = 0.2; lS = 0.05;

// Cells x-direction
nbxLE    = 5;
nbxLConv = 20;
nbxLCol  = 5;
nbxLDiv  = 20;
nbxLS    = 5;

// Cells y-direction
nby = 10;

// Boundary condition numbers
bcAxis = 1;
bcWall = 2;
bcInlet = 3;
bcOutlet = 4;
surfFluid = 10;

//----------------- GEO -----------------

// Geometry inlet
Point(1) = {0, 0, 0};
Point(2) = {lE, 0, 0};
Point(3) = {lE, rE, 0};
Point(4) = {0, rE, 0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Geometry convergent
Point(5) = {lE+lConv, 0, 0};
Point(6) = {lE+lConv, rC, 0};
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,3};
Line Loop(2) = {-2,5,6,7};
Plane Surface(2) = {2};

// Geometry throat
Point(7) = {lE+lConv+lCol, 0, 0};
Point(8) = {lE+lConv+lCol, rC, 0};
Line(8) = {5,7};
Line(9) = {7,8};
Line(10) = {8,6};
Line Loop(3) = {-6,8,9,10};
Plane Surface(3) = {3};

// Geometry divergent
Point(9) = {lE+lConv+lCol+lDiv, 0, 0};
Point(10) = {lE+lConv+lCol+lDiv, rS, 0};
Line(11) = {7,9};
Line(12) = {9,10};
Line(13) = {10,8};
Line Loop(4) = {-9,11,12,13};
Plane Surface(4) = {4};

// Geometry outlet
Point(11) = {lE+lConv+lCol+lDiv+lS, 0, 0};
Point(12) = {lE+lConv+lCol+lDiv+lS, rS, 0};
Line(14) = {9,11};
Line(15) = {11,12};
Line(16) = {12,10};
Line Loop(5) = {-12,14,15,16};
Plane Surface(5) = {5};

// Structured mesh
Transfinite Line {1,3} = nbxLE;
Transfinite Line {5,7} = nbxLConv;
Transfinite Line {8,10} = nbxLCol;
Transfinite Line {11,13} = nbxLDiv;
Transfinite Line {14,16} = nbxLS;
Transfinite Line {2,4,6,9,12,15} = nby;
Transfinite Surface "*";
Recombine Surface "*";

// Boundary conditions and fluid domain
Physical Line(bcAxis) = {1,5,8,11,14};
Physical Line(bcWall) = {3,7,10,13,16};
Physical Line(bcInlet) = {4};
Physical Line(bcOutlet) = {15};
Physical Surface(surfFluid) = {1,2,3,4,5};
