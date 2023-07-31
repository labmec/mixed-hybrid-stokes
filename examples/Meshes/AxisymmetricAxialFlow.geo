elr = 5;
elz = 5;
ndivr = elr + 1;
ndivz = elz + 1;

Ri = 2.0;
Re = 4.0;
z = -1.0;
L = 2.0;

Point(1) = {Ri, z, 0.0, 1.0};
//+
Point(2) = {Re, z, 0.0, 1.0};
//+
Point(3) = {Re, z+L, 0.0, 1.0};
//+
Point(4) = {Ri, z+L, 0.0, 1.0};
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

Physical Curve("NormalStress", 2) = {1,3};
//+
Physical Curve("NormalVelocity", 3) = {2,4};
//+
Physical Curve("OuterTangentialVelocity", 4) = {2};
//+
Physical Curve("InnerTangentialVelocity", 5) = {4};
//+
Physical Curve("NoTangentialVelocity", 6) = {1,3};
//+

Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Curve {4} = ndivz Using Progression 1;
//+
Transfinite Curve {1} = ndivr Using Progression 1;
//+
Transfinite Curve {2} = ndivz Using Progression 1;
//+
Transfinite Curve {3} = ndivr Using Progression 1;
//+

Recombine Surface {1};
