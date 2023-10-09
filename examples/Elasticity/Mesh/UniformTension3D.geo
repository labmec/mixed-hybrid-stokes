// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

el = 2;
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


//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, -L} {
  Surface{1}; 
}

Physical Volume("Domain", 1) = {1};
//+
Physical Surface("UnitNormalStress", 2) = {5};
//+
Physical Surface("ZeroNormalStress", 3) = {2};
//+
Physical Surface("ZeroNormalDisplacement", 4) = {6, 3, 4, 1};
//+
Physical Surface("ZeroTangentialStress", 5) = {3, 2, 6, 5, 1, 4};
//+

Transfinite Curve {4, 6, 7, 5, 3, 9, 12, 2, 1, 8, 11, 10} = ndiv Using Progression 1;
//+
Transfinite Surface {2} = {6, 4, 3, 5};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Surface {5} = {2, 8, 5, 3};
//+
Transfinite Surface {3} = {1, 4, 6, 7};
//+
Transfinite Surface {4} = {1, 7, 8, 2};
//+
Transfinite Surface {6} = {6, 5, 8, 7};
//+
Recombine Surface {2, 3, 6, 5, 1, 4};
//+

Transfinite Volume{1} = {1, 2, 3, 4, 7, 8, 5, 6};
