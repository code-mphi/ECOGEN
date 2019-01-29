// lc = 1. ;
// xx = 0;

// // Splines (CatmullRom)
// p = newp;
// Point(p) = {xx,  0,  0, 0.3*lc} ;
// Point(p+1) = {xx+.5, 0,  0, lc} ;
// Point(p+2) = {xx+.5, 1, 0, 0.1*lc} ;
// Point(p+3) = {xx, 1, 0, lc} ;
// l = newreg;
// Spline(l) = {p+3,p+2,p+1,p};
// Line(l+1) = {p,p+3};
// s = newreg;
// Line Loop(s) = {-l,-(l+1)};
// Plane Surface(s+1) = s;

// // Duplicate the surfaces, and use uniform mesh
// p1 = newp;
// Translate {0,-1.5,0} {
//   Duplicata { Surface{6:18:4}; }
// }
// p2 = newp;
// Printf("p1 p2 = %g %g", p1, p2);

// Characteristic Length {p1:p2-1} = lc/5 ;


// Gmsh project created on Fri Sep 12 15:39:15 2014
k=0.2;
Point(1) = {0, 0, 0, k};
Point(2) = {2, 0, 0, k};
Point(3) = {3, 2, 0, k};
Point(4) = {1, 2, 0, k};
Point(5) = {1.5, 1, 0, k};

Spline(2) = {4, 5, 1, 2};
Spline(3) = {2, 3, 4};
Line Loop(3) = {2,3};
Plane Surface(1) = {3};

Extrude {0, 0, 1} {
  Surface{1};
}
Physical Surface(1) = {14, 15};
Physical Surface(2) = {1, 10};
Physical Volume(17) = {1};

Recombine Surface {14, 15, 10, 1};
