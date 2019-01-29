// Maillage Tuyere 2D simple
//
//    rE _                                          
//    rS _ |\__  _  _  _  _  _   
//         |   \__          ___/|
//         |      \__   ___/    |
//    rC _ | _  _  _ \ /        |
//     0 _ |__________|_________|
//
//         <---------><-------->
//             lConv     lDiv    
//
//

//Maillage Tuyere 2D simple - Details geometrie

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

//---------------------------------DONNEES------------------------------------
//Dimensions principales
rE = 0.14657; rC = 0.06406; rS = 0.14657;
lConv = 0.5; lDiv = 0.5;
//Number de mailles selon X
nbxLConv = 50;
nbxLDiv  = 50;
//number de mailles selon Y
nby = 1;
//Numero des conditions aux limites
condLimAxe    = 1;
condLimParoi  = 2;
condLimEntree = 3;
condLimSortie = 4;
surfaceFluide = 10;
//----------------------------------------------------------------------------

//Geometrie Entree
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
//Geometrie Divergent
Point(5) = {lConv+lDiv, 0, 0};
Point(6) = {lConv+lDiv, rS, 0};
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,3};
Line Loop(2) = {5,6,7,-2};
Plane Surface(2) = {2};
//Maillages
Transfinite Line {1,3} = nbxLConv;
Transfinite Line {5,7} = nbxLDiv;
Transfinite Line {2,4,6} = nby;
Transfinite Surface "*";
Recombine Surface "*";
//Domaines Physiques et CL
Physical Line(condLimAxe) = {1,5};
Physical Line(condLimParoi) = {3,7};
Physical Line(condLimEntree) = {4};
Physical Line(condLimSortie) = {6};
Physical Surface(surfaceFluide) = {1,2};
