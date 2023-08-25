elr = 45;
elz = 50;
ndivr = elr + 1;
ndivz = elz + 1;
ndivo = 2*elr/15 + 1;

Ri = 0.001;
Re = 1.5;
z = 0.0;
Lin = 5.0;
Lout = 5.0;

Point(1) = {Ri, -Lin, 0.0, 1.0};
//+
Point(2) = {Re, -Lin, 0.0, 1.0};
//+
Point(3) = {Re, Lout, 0.0, 1.0};
//+
Point(4) = {Ri, Lout, 0.0, 1.0};
//+
Point(5) = {Ri, 0.0, 0.0, 1.0};
//+
Point(6) = {1*Re/15, 0.0, 0.0, 1.0};
//+
Point(7) = {3*Re/15, 0.0, 0.0, 1.0};
//+
Point(8) = {5*Re/15, 0.0, 0.0, 1.0};
//+
Point(9) = {7*Re/15, 0.0, 0.0, 1.0};
//+
Point(10) = {9*Re/15, 0.0, 0.0, 1.0};
//+
Point(11) = {11*Re/15, 0.0, 0.0, 1.0};
//+
Point(12) = {13*Re/15, 0.0, 0.0, 1.0};
//+
Point(13) = {Re, 0.0, 0.0, 1.0};
//+
Point(14) = {Re, 0.0, 0.0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 13};
//+
Line(3) = {14, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 1};
//+
Line(7) = {13, 12};
//+
Line(8) = {12, 11};
//+
Line(9) = {11, 10};
//+
Line(10) = {10, 9};
//+
Line(11) = {9, 8};
//+
Line(12) = {8, 7};
//+
Line(13) = {7, 6};
//+
Line(14) = {6, 5};
//+
Line(15) = {6, 7};
//+
Line(16) = {8, 9};
//+
Line(17) = {10, 11};
//+
Line(18) = {12, 14};
//+

Curve Loop(1) = {1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, -14, 15, -12, 16, -10, 17, -8, 18};
//+
Plane Surface(2) = {2};
//+

Physical Surface("Domain", 1) = {1,2};

Physical Curve("InFlow", 2) = {1};
//+
Physical Curve("OutFlow", 3) = {4};
//+
Physical Curve("NoNormalFlow", 4) = {2,3,7,9,11,13,15,16,17,18};
//+
Physical Curve("NoSlipWall", 5) = {1,2,3,4,7,9,11,13,15,16,17,18};
//+
Physical Curve("Infinitesimal_Tube", 6) = {5,6};
//+
Physical Point("InflowTube", 7) = {1};
//+
Physical Point("OutflowTube", 8) = {4};
//+

//Transfinite Surface {1} = {1, 2, 13, 5};
//+
//Transfinite Surface {2} = {3, 4, 5, 14};
//+
Transfinite Curve {1} = 5 Using Progression 1;
//+
Transfinite Curve {2} = 51 Using Progression 0.909;
//+
Transfinite Curve {3} = 51 Using Progression 1.1;
//+
Transfinite Curve {4} = 5 Using Progression 1;
//+
Transfinite Curve {5} = 51 Using Progression 0.909;
//+
Transfinite Curve {6} = 51 Using Progression 1.1;
//+
Transfinite Curve {7} = 21 Using Progression 1;
//+
Transfinite Curve {8} = 21 Using Progression 1;
//+
Transfinite Curve {9} = 21 Using Progression 1;
//+
Transfinite Curve {10} = 21 Using Progression 1;
//+
Transfinite Curve {11} = 21 Using Progression 1;
//+
Transfinite Curve {12} = 21 Using Progression 1;
//+
Transfinite Curve {13} = 21 Using Progression 1;
//+
Transfinite Curve {14} = 11 Using Progression 1;
//+
Transfinite Curve {15} = 21 Using Progression 1;
//+
Transfinite Curve {16} = 21 Using Progression 1;
//+
Transfinite Curve {17} = 21 Using Progression 1;
//+
Transfinite Curve {18} = 21 Using Progression 1;
//+

Recombine Surface {1};
Recombine Surface {2};
//+
