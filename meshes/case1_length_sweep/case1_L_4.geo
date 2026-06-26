// problem_domain_2D.geo
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// --------------------
// Parameters
// --------------------
Lb = 4.6;           // beam length [m]
Ld = 30.0;          // damping zone length [m]
Lf = 60.0;          // central fluid region length [m]
L  = 2*Ld + Lf;     // total tank length [m] = 75.0
h  = 19.0;          // water depth [m]
lc_fs   = 0.15;     // free surface and structure [m]
lc_damp = 1.00;     // damping zones [m]
lc_deep = 2.00;     // deep region [m]
lc_far  = 2.00;     // bottom and side walls [m]
z_deep  = -5.0;     // depth below which mesh is coarsened [m]

// Damping zone boundaries
xdin  = Ld;         // 20.0
xdout = Ld + Lf;    // 55.0

// Beam coordinates (centered in the middle region)
xL = Ld + (Lf - Lb)/2;   // 31.25
xR = xL + Lb;             // 43.75

// Coordinates in (x,z) with free surface at z=0 and bottom at z=-h
Point(1) = {0,    -h, 0, lc_far};
Point(2) = {L,    -h, 0, lc_far};
Point(3) = {L,     0, 0, lc_damp};
Point(4) = {0,     0, 0, lc_damp};
Point(5) = {xL,    0, 0, lc_fs};
Point(6) = {xR,    0, 0, lc_fs};
Point(9) = {xdin,  0, 0, lc_fs};
Point(10)= {xdout, 0, 0, lc_fs};

// --------------------
// Boundary curves (counterclockwise)
// --------------------
Line(1)  = {1,2};    // Bottom
Line(2)  = {2,3};    // Outlet

Line(11) = {3,10};   // FreeSurface outlet damping
Line(3)  = {10,6};   // FreeSurface right of beam
Line(5)  = {6,5};    // Structure
Line(6)  = {5,9};    // FreeSurface left of beam
Line(12) = {9,4};    // FreeSurface inlet damping

Line(4)  = {4,1};    // Inlet

// --------------------
// Surface
// --------------------
Curve Loop(10) = {1,2,11,3,5,6,12,4};
Plane Surface(11) = {10};

// --------------------
// Refinement fields
// --------------------
Field[1] = Box;
Field[1].XMin = 0;
Field[1].XMax = xdin;
Field[1].YMin = -h;
Field[1].YMax = 0;
Field[1].VIn  = lc_damp;
Field[1].VOut = lc_fs;

Field[2] = Box;
Field[2].XMin = xdout;
Field[2].XMax = L;
Field[2].YMin = -h;
Field[2].YMax = 0;
Field[2].VIn  = lc_damp;
Field[2].VOut = lc_fs;

Field[3] = Box;
Field[3].XMin = 0;
Field[3].XMax = L;
Field[3].YMin = -h;
Field[3].YMax = z_deep;
Field[3].VIn  = lc_deep;
Field[3].VOut = lc_fs;

Field[4] = Max;
Field[4].FieldsList = {1, 2, 3};

Background Field = 4;

// --------------------
// Physical groups (TAGS)
// --------------------
Physical Curve("Bottom")      = {1};
Physical Curve("Outlet")      = {2};
Physical Curve("FreeSurface") = {11,3,6,12};
Physical Curve("Inlet")       = {4};
Physical Curve("Structure")   = {5};

Physical Surface("Fluid")     = {11};

Physical Point("BottomLeft") = {1};
Physical Point("DampingInlet")  = {9};
Physical Point("DampingOutlet") = {10};
Physical Point("BeamLeft")      = {5};
Physical Point("BeamRight")     = {6};
