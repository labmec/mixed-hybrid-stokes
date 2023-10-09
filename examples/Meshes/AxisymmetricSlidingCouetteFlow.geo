elr = 16;
elz = 16;
ndivr = elr + 1;
ndivz = elz + 1;
//+

R1 = 1.0;
R2 = 2.0;
z = 0.0;
L = 1.0;
//+

Point(1) = {R1, z, 0, 1.0};
//+
Point(2) = {R2, z, 0, 1.0};
//+
Point(3) = {R2, z+L, 0, 1.0};
//+
Point(4) = {R1, z+L, 0, 1.0};
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

Physical Curve("UnitNormalStress", 2) = {1,3};
//+
Physical Curve("NoNormVelocity", 3) = {2,4};
//+
Physical Curve("NoTanVelocity", 4) = {1,2,3};
//+
Physical Curve("TanVelocity", 5) = {4};

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
