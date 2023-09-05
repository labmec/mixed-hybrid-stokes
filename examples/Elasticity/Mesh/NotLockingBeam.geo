// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

elx = 7;
ely = 2;
ndivx = elx + 1;
ndivy = ely + 1;

L = 10;
h = 1;

Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {L, 0.0, 0.0, 1.0};
//+
Point(3) = {L, h, 0.0, 1.0};
//+
Point(4) = {0, h, 0.0, 1.0};
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

Physical Curve("ZeroNormalVelocity", 2) = {4};
//+
Physical Curve("ZeroNormalStress", 3) = {1,2,3};
//+
Physical Curve("ZeroTangentialVelocity", 4) = {4};
//+
Physical Curve("ZeroTangentialStress", 5) = {1,3};
//+
Physical Curve("TangentialStress", 6) = {2};
//+

Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Curve {4} = ndivy Using Progression 1;
//+
Transfinite Curve {1} = ndivx Using Progression 1;
//+
Transfinite Curve {2} = ndivy Using Progression 1;
//+
Transfinite Curve {3} = ndivx Using Progression 1;
//+

Recombine Surface {1};
