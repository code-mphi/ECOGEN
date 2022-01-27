// NACA 66-312mod geometry

// Default unit = m

// --- Geometry parameters ---
aoa = 6.;        // Angle of attack (°) (> 0 to get leading edge moved up)
chord = 150.e-3; // Length from leading edge to trailing edge
lx = 1.95;       // Lenght of the domain
ly = 192.e-3;    // Height of the domain

// dx1 = 0.0025; // 3k4
// dx2 = 0.025;

dx1 = 0.0015; // 8k9
dx2 = 0.015;

// dx1 = 0.00125; // 12k
// dx2 = 0.0125;

// dx1 = 0.001; // 19k
// dx2 = 0.01;

// --- Boundary conditions ---
boundCondInflow = 1;
boundCondOutflow = 2;
boundCondWalls = 3;
boundCondFoil = 4;
fluidSurface = 10;

// --- Geometry ---

// Set angle in the flow direction and convert to radian
aoa = - aoa * Pi / 180; 

lx_lead_edge = 0.3 * lx;
ly_lead_edge = 0.5 * ly;

// Create surrouding box
Point(newp) = {0., 0., 0., dx2};
Point(newp) = {0., ly, 0., dx2};
Point(newp) = {lx, ly, 0., dx2};
Point(newp) = {lx, 0., 0., dx2};

For i In {1:3}
  Line(newl) = {i, i+1};
EndFor
Line(newl) = {4, 1};

// Extrados
x_pos_extra = {
-0.000009,
-0.000124,
0.000401,
0.001576,
0.003445,
0.006158,
0.010008,
0.015387,
0.022659,
0.032046,
0.043639,
0.057431,
0.073337,
0.09124,
0.111019,
0.132551,
0.155704,
0.180355,
0.206385,
0.233668,
0.262076,
0.29148,
0.321756,
0.352777,
0.384417,
0.41655,
0.449049,
0.481787,
0.514635,
0.547462,
0.58014,
0.612539,
0.644532,
0.675988,
0.706782,
0.736788,
0.76588,
0.793924,
0.820781,
0.846343,
0.870518,
0.893207,
0.914308,
0.933718,
0.951337,
0.967069,
0.980845,
0.992676,
1.002725,
1.011376,
1.01924
};

y_pos_extra = {
0.000051,
0.002603,
0.00515,
0.007641,
0.010119,
0.012727,
0.015633,
0.018917,
0.022569,
0.026559,
0.030809,
0.035195,
0.039634,
0.044059,
0.048395,
0.05259,
0.05663,
0.060466,
0.064017,
0.067244,
0.070148,
0.072722,
0.074937,
0.076761,
0.078178,
0.079199,
0.079756,
0.079768,
0.079245,
0.078206,
0.076647,
0.074584,
0.07203,
0.069001,
0.065526,
0.061649,
0.057409,
0.052776,
0.047734,
0.042421,
0.03703,
0.031692,
0.026511,
0.021589,
0.017018,
0.012873,
0.009205,
0.00603,
0.003316,
0.000969,
-0.001176
};

p_extra_b = newp;

For i In {0:#x_pos_extra[]-1}
  Point(p_extra_b + i) = {Cos(aoa) * x_pos_extra[i] * chord - Sin(aoa) * y_pos_extra[i] * chord + lx_lead_edge,
    Sin(aoa) * x_pos_extra[i] * chord + Cos(aoa) * y_pos_extra[i] * chord + ly_lead_edge, 
    0.,
    dx1
  };
  p_extra_e = p_extra_b + i;
EndFor

// Intrados
x_pos_intra = {
1.011176,
1.002356,
0.992151,
0.980162,
0.966219,
0.950309,
0.932506,
0.912916,
0.891656,
0.868841,
0.844591,
0.819029,
0.792279,
0.764464,
0.735692,
0.706072,
0.675724,
0.644765,
0.613311,
0.581481,
0.549393,
0.517169,
0.484928,
0.45279,
0.420875,
0.389304,
0.358197,
0.327672,
0.29785,
0.268849,
0.240789,
0.213789,
0.187971,
0.163454,
0.140356,
0.118796,
0.098896,
0.080778,
0.064561,
0.05036,
0.038265,
0.028302,
0.020396,
0.014349,
0.009834,
0.006473,
0.003946,
0.002056,
0.000731
};

y_pos_intra = {
-0.002367,
-0.003627,
-0.005026,
-0.006587,
-0.008294,
-0.010109,
-0.011999,
-0.013934,
-0.015881,
-0.017807,
-0.019709,
-0.021624,
-0.02365,
-0.02583,
-0.02803,
-0.030126,
-0.032114,
-0.033964,
-0.03561,
-0.037033,
-0.038221,
-0.039155,
-0.039828,
-0.04022,
-0.040311,
-0.040166,
-0.039849,
-0.039343,
-0.038643,
-0.037773,
-0.036747,
-0.035556,
-0.034178,
-0.032624,
-0.030944,
-0.029164,
-0.027272,
-0.025273,
-0.023196,
-0.021057,
-0.018879,
-0.016722,
-0.014623,
-0.012586,
-0.010627,
-0.008721,
-0.006775,
-0.004687,
-0.002407
};

p_intra_b = newp;

For i In {0:#x_pos_intra[]-1}
  Point(p_intra_b + i) = {Cos(aoa) * x_pos_intra[i] * chord - Sin(aoa) * y_pos_intra[i] * chord + lx_lead_edge,
    Sin(aoa) * x_pos_intra[i] * chord + Cos(aoa) * y_pos_intra[i] * chord + ly_lead_edge, 
    0.,
    dx1
  };
  p_intra_e = p_intra_b + i;
EndFor

// Create foil shape
l_foil_b = newl;
For i In {p_extra_b:p_intra_e-1}
  Line(l_foil_b + i-p_extra_b) = {i, i+1};
  l_foil_e = l_foil_b + i-p_extra_b;
EndFor
Line(newl) = {p_intra_e, p_extra_b};

// Build fluid surface
Line Loop(1) = {1:4}; // Loop domain
Line Loop(2) = {l_foil_b:l_foil_e+1}; // Loop foil
Plane Surface(1) = {1, 2};

// Set boundaries and volume
Physical Line(boundCondInflow) = {1};
Physical Line(boundCondOutflow) = {3};
Physical Line(boundCondWalls) = {2,4};
Physical Line(boundCondFoil) = {l_foil_b:l_foil_e+1};
Physical Surface(fluidSurface) = {1};
