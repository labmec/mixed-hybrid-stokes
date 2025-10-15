// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

el = 16;
ndiv = el + 1;


Point(1) = {-0.5, 0.0, 0.0, 1.0};
//+
Point(2) = {1.5, 0.0, 0.0, 1.0};
//+
Point(3) = {1.5, 2.0, 0.0, 1.0};
//+
Point(4) = {-0.5, 2.0, 0.0, 1.0};
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

Physical Curve("BottomNormal", 2) = {1};
//+
Physical Curve("RightNormal", 3) = {2};
//+
Physical Curve("TopNormal", 4) = {3};
//+
Physical Curve("LeftNormal", 5) = {4};
//+

Physical Curve("BottomTan", 6) = {1};
//+
Physical Curve("RightTan", 7) = {2};
//+
Physical Curve("TopTan", 8) = {3};
//+
Physical Curve("LeftTan", 9) = {4};
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
