elr = 20;
elz = 20;
ndivr = elr + 1;
ndivz = elz + 1;
ndivo = 0.5*elr + 1;

Ri = 0.001;
Re = 1.0;
z = 0.0;
L = 5.0;
OR = 0.5*Re;

Point(1) = {Ri, z, 0.0, 1.0};
//+
Point(2) = {Re, z, 0.0, 1.0};
//+
Point(3) = {Re, z+(L/2), 0.0, 1.0};
//+
Point(4) = {Re-OR, z+(L/2), 0.0, 1.0};
//+
Point(5) = {Re, z+(L/2), 0.0, 1.0};
//+
Point(6) = {Re, z+L, 0.0, 1.0};
//+
Point(7) = {Ri, z+L, 0.0, 1.0};
//+
Point(8) = {Ri, z+(L/2), 0.0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 1};
//+
Line(9) = {4, 8};
//+

Curve Loop(1) = {1, 2, 3, 9, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, -9, 4};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Domain", 1) = {1,2};

Physical Curve("InFlow", 2) = {1};
//+
Physical Curve("OutFlow", 3) = {6};
//+
Physical Curve("NoNormalFlow", 4) = {2,3,4,5};
//+
Physical Curve("NoSlipWall", 5) = {1,2,3,4,5,6};
//+
Physical Curve("Infinitesimal_Tube", 6) = {7,8};
//+
Physical Point("InflowTube", 7) = {1};
//+
Physical Point("OutflowTube", 8) = {7};
//+

Transfinite Surface {1} = {1, 2, 3, 8};
//+
Transfinite Surface {2} = {8, 5, 6, 7};
//+
Transfinite Curve {1} = ndivr Using Progression 1;
//+
Transfinite Curve {2} = ndivz Using Progression 1;
//+
Transfinite Curve {3} = ndivo Using Progression 1;
//+
Transfinite Curve {4} = ndivo Using Progression 1;
//+
Transfinite Curve {5} = ndivz Using Progression 1;
//+
Transfinite Curve {6} = ndivr Using Progression 1;
//+
Transfinite Curve {7} = ndivz Using Progression 1;
//+
Transfinite Curve {8} = ndivz Using Progression 1;
//+
Transfinite Curve {9} = ndivr-ndivo+1 Using Progression 1;
//+

Recombine Surface {1};
Recombine Surface {2};
