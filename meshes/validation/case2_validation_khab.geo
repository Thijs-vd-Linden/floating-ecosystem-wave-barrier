SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// --------------------
// Parameters
// --------------------
Lb = 12.5;
Ld = 50.0;
Lf = 50.0;
L  = 2*Ld + Lf;
H  = 1.1;
dp = 0.75;

lc = 0.20;

xdin  = Ld;
xdout = Ld + Lf;

xL = Ld + (Lf - Lb)/2;
xR = xL + Lb;

// --------------------
// Points
// --------------------
Point(1) = {0,   0, 0, lc};
Point(2) = {L,   0, 0, lc};
Point(3) = {L,   H, 0, lc};
Point(4) = {0,   H, 0, lc};

Point(9)  = {xdin,  H, 0, lc};
Point(10) = {xdout, H, 0, lc};

Point(5) = {xL,  H,    0, lc};
Point(6) = {xR,  H,    0, lc};

Point(7) = {xL,  H-dp, 0, lc};
Point(8) = {xR,  H-dp, 0, lc};

// --------------------
// Boundary curves
// --------------------
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {4,1};

Line(4)  = {4,  9};
Line(10) = {9,  5};
Line(5)  = {6, 10};
Line(11) = {10, 3};

Line(12) = {5,6};

Line(7)  = {5,7};
Line(8)  = {7,8};
Line(9)  = {8,6};

// --------------------
// Surfaces
// --------------------
Curve Loop(10) = {1, 2, -11, -5, -9, -8, -7, -10, -4, 3};
Plane Surface(11) = {10};

Curve Loop(20) = {-12, 7, 8, 9};
Plane Surface(21) = {20};

// --------------------
// Physical groups
// --------------------
Physical Curve("Bottom")      = {1};
Physical Curve("Outlet")      = {2};
Physical Curve("Inlet")       = {3};
Physical Curve("FreeSurface") = {4, 10, 5, 11};
Physical Curve("Structure")   = {12};

Physical Surface("Fluid")     = {11};
Physical Surface("Porous")    = {21};

Physical Point("BottomLeft")       = {1};
Physical Point("FreeSurfaceLeft")  = {4};
Physical Point("FreeSurfaceRight") = {3};
Physical Point("DampingInlet")     = {9};
Physical Point("DampingOutlet")    = {10};
Physical Point("BeamLeft")         = {5};
Physical Point("BeamRight")        = {6};