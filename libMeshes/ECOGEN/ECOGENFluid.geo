dx = 1.;

//Whole domain
Point(1) = {-3, -3, 0, dx};
Point(4) = {20, -3, 0, dx};
Point(5) = {20, 8, 0, dx};
Point(6) = {-3, 8, 0, dx};
Line(1) = {1, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};
Line Loop(1) = {3, 4, 1, 2};

//Letter E
Point(7) = {0.75, 1, 0, dx};
Point(8) = {0.75, 4, 0, dx};
Point(9) = {2.5, 4, 0, dx};
Point(10) = {2.5, 1, 0, dx};
Point(11) = {1.25, 1.5, 0, dx};
Point(12) = {2.5, 1.5, 0, dx};
Point(13) = {2.5, 3.5, 0, dx};
Point(14) = {1.25, 3.5, 0, dx};
Point(15) = {1.25, 2.25, 0, dx};
Point(16) = {1.25, 2.75, 0, dx};
Point(17) = {2., 2.75, 0, dx};
Point(18) = {2., 2.25, 0, dx};
Line(5) = {7, 10};
Line(6) = {10, 12};
Line(7) = {12, 11};
Line(8) = {11, 15};
Line(9) = {15, 18};
Line(10) = {18, 17};
Line(11) = {17, 16};
Line(12) = {16, 14};
Line(13) = {14, 13};
Line(14) = {13, 9};
Line(15) = {9, 8};
Line(16) = {8, 7};
Line Loop(2) = {15, 16, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

//Letter C
Point(19) = {4.5, 1, 0, dx};
Point(20) = {4.5, 1.5, 0, dx};
Point(21) = {4.5, 4, 0, dx};
Point(22) = {4.5, 3.5, 0, dx};
Point(23) = {4.5, 2.5, 0, dx};
Circle(20) = {21, 23, 19};
Circle(21) = {22, 23, 20};
Point(24) = {4.75, 1.5, 0, dx};
Point(25) = {4.75, 1., 0, dx};
Point(26) = {4.75, 3.5, 0, dx};
Point(27) = {4.75, 4, 0, dx};
Line(22) = {21, 27};
Line(23) = {27, 26};
Line(24) = {26, 22};
Line(25) = {20, 24};
Line(26) = {24, 25};
Line(27) = {25, 19};
Line Loop(3) = {22, 23, 24, 21, 25, 26, 27, -20};

//Letter O
//Ellipse ( expression ) = { expression, expression, expression, expression <, ...> };
//Creates an ellipse arc. The four expressions on the right-hand-side define the start point, the center point, a major axis point and the end point of the ellipse.
Point(28) = {6.5, 1, 0, dx};
Point(29) = {6.5, 2.5, 0, dx};
Point(30) = {6.5, 4, 0, dx};
Point(31) = {5.25, 2.5, 0, dx};
Point(32) = {7.75, 2.5, 0, dx};
Ellipse(28) = {28, 29, 30, 31};
Ellipse(29) = {31, 29, 32, 30};
Ellipse(30) = {30, 29, 28, 32};
Ellipse(31) = {32, 29, 31, 28};
Point(33) = {6.5, 1.5, 0, dx};
Point(34) = {6.5, 3.5, 0, dx};
Point(35) = {5.75, 2.5, 0, dx};
Point(36) = {7.25, 2.5, 0, dx};
Ellipse(32) = {33, 29, 34, 35};
Ellipse(33) = {35, 29, 36, 34};
Ellipse(34) = {34, 29, 33, 36};
Ellipse(35) = {36, 29, 35, 33};
Line Loop(4) = {29, 30, 31, 28};
Line Loop(5) = {33, 34, 35, 32};

//Letter G
Point(37) = {8.25, 2.5, 0, dx};
Point(38) = {9.75, 2.5, 0, dx};
Point(39) = {9.75, 1, 0, dx};
Point(40) = {9.75, 4, 0, dx};
Point(41) = {9.75, 1.5, 0, dx};
Point(42) = {9.75, 3.5, 0, dx};
Point(43) = {8.75, 2.5, 0, dx};
Circle(36) = {40, 38, 39};
Circle(37) = {42, 38, 41};
Point(44) = {10.25, 4, 0, dx};
Point(45) = {10.25, 3.5, 0, dx};
Line(38) = {44, 45};
Line(39) = {44, 40};
Line(40) = {42, 45};
Point(46) = {10.25, 1, 0, dx};
Point(47) = {10.25, 2.5, 0, dx};
Point(48) = {9.75, 2., 0, dx};
Point(49) = {9.5, 2., 0, dx};
Point(50) = {9.5, 2.5, 0, dx};
Line(41) = {46, 39};
Line(42) = {47, 46};
Line(43) = {47, 50};
Line(44) = {50, 49};
Line(45) = {49, 48};
Line(46) = {48, 41};
Line Loop(6) = {36, -41, -42, 43, 44, 45, 46, -37, 40, -38, 39};

//Letter E
x = 10;
Point(51) = {0.75+x, 1, 0, dx};
Point(52) = {0.75+x, 4, 0, dx};
Point(53) = {2.5+x, 4, 0, dx};
Point(54) = {2.5+x, 1, 0, dx};
Point(55) = {1.25+x, 1.5, 0, dx};
Point(56) = {2.5+x, 1.5, 0, dx};
Point(57) = {2.5+x, 3.5, 0, dx};
Point(58) = {1.25+x, 3.5, 0, dx};
Point(59) = {1.25+x, 2.25, 0, dx};
Point(60) = {1.25+x, 2.75, 0, dx};
Point(61) = {2.+x, 2.75, 0, dx};
Point(62) = {2.+x, 2.25, 0, dx};
Line(59) = {51, 54};
Line(60) = {54, 56};
Line(61) = {56, 55};
Line(62) = {55, 59};
Line(63) = {59, 62};
Line(64) = {62, 61};
Line(65) = {61, 60};
Line(66) = {60, 58};
Line(67) = {58, 57};
Line(68) = {57, 53};
Line(69) = {53, 52};
Line(70) = {52, 51};
Line Loop(7) = {69, 70, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68};

//Letter N
Point(63) = {13, 1, 0, dx};
Point(64) = {13, 4, 0, dx};
Point(65) = {13.5, 1, 0, dx};
Point(66) = {13.5, 4, 0, dx};
Point(67) = {13.5, 3., 0, dx};
Point(68) = {15, 1, 0, dx};
Point(69) = {15, 4, 0, dx};
Point(70) = {14.5, 1, 0, dx};
Point(71) = {14.5, 4, 0, dx};
Point(72) = {14.5, 2, 0, dx};
Line(71) = {63, 64};
Line(72) = {64, 66};
Line(73) = {66, 72};
Line(74) = {72, 71};
Line(75) = {71, 69};
Line(76) = {69, 68};
Line(77) = {68, 70};
Line(78) = {70, 67};
Line(79) = {67, 65};
Line(80) = {65, 63};
Line Loop(8) = {72, 73, 74, 75, 76, 77, 78, 79, 80, 71};

//Physical lines and surface
Physical Line(1) = {1,2,4}; //Walls
Physical Line(2) = {3}; //output

Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8}; //Domaine Gazeux
Plane Surface(2) = {5}; //interieur du O
Physical Surface(10) = {1,2};
Plane Surface(3) = {2}; 
Plane Surface(4) = {3};
Plane Surface(5) = {4,5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

Physical Surface(11) = {3}; //E
Physical Surface(12) = {4}; //C
Physical Surface(13) = {5}; //O
Physical Surface(14) = {6}; //G
Physical Surface(15) = {7}; //E
Physical Surface(16) = {8}; //N

//Physical Surface(11) = {3,4,5,6,7,8};



