// Geometry Venturi 8°
// Data from Le Martelot et al. 2013

// Geometry parameters (m)
xA = 0.;    yA = 0.;
xB = 0.1;   yB = 0.;
xC = 0.153; yC = 0.0157;
xD = 0.588; yD = -0.0517;
xE = 1.225; yE = -0.114;
xF = 0.;    yF = 0.0488;
xG = 0.271; yG = 0.0488;
xH = 1.233; yH = -0.00845;

// Mesh parameters
// 12k 
// dsBig   = 0.008;
// dsMid   = 0.003;
// dsSmall = 0.003; 
// 21k 
// dsBig   = 0.005;
// dsMid   = 0.0025;
// dsSmall = 0.0025; 
// 42k 
dsBig   = 0.005;
dsMid   = 0.005;
dsSmall = 0.001; 
// 95k
// dsBig   = 0.0035;
// dsMid   = 0.0035;
// dsSmall = 0.00065; 
// 167k
// dsBig   = 0.0025;
// dsMid   = 0.0025;
// dsSmall = 0.0005; 
// 387k
// dsBig   = 0.002;
// dsMid   = 0.002;
// dsSmall = 0.0003; 

// Boundary conditions
boundCondInflow = 1;
boundCondOutflow = 2;
boundCondWalls = 3;

// Geometry
// Be aware that points ordering follows
// the convention used in the paper
Point(1) = {xA, yA, 0, dsSmall};
Point(2) = {xB, yB, 0, dsSmall};
Point(3) = {xC, yC, 0, dsSmall};
Point(4) = {xD, yD, 0, dsMid};
Point(5) = {xE, yE, 0, dsBig};
Point(6) = {xF, yF, 0, dsSmall};
Point(7) = {xG, yG, 0, dsSmall};
Point(8) = {xH, yH, 0, dsBig};

// Line definition : A-B-C-D-E-H-G-F-A
For i In {1:4}
    Line(newl) = {i, i+1};
EndFor
Line(newl) = {5, 8};
Line(newl) = {8, 7};
Line(newl) = {7, 6};
Line(newl) = {6, 1}; 

Line Loop(1) = {1:8};
Plane Surface(1) = {1};

// Boundary conditions and physical surface
Physical Line(boundCondInflow) = {8};
Physical Line(boundCondOutflow) = {5};
Physical Line(boundCondWalls) = {1:4, 6, 7};
Physical Surface(10) = {1};

