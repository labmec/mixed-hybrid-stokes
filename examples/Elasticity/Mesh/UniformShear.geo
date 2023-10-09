// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

el = 1;
ndiv = el + 1;

L = 1.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {L, 0.0, 0.0, 1.0};
//+
Point(3) = {L, L, 0.0, 1.0};
//+
Point(4) = {0, L, 0.0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Domain", 1) = {1};

Physical Curve("ZeroNormalDisplacement", 2) = {1};
//+
Physical Curve("ZeroNormalStress", 3) = {2,3,4};
//+
Physical Curve("ZeroTangentialDisplacement", 4) = {1,2,4};
//+
Physical Curve("UnitTangentialDisplacement", 5) = {3};
//+

Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Curve {4} = ndiv Using Progression 1;
//+
Transfinite Curve {1} = ndiv Using Progression 1;
//+
Transfinite Curve {2} = ndiv Using Progression 1;
//+
Transfinite Curve {3} = ndiv Using Progression 1;
//+

Recombine Surface {1};
