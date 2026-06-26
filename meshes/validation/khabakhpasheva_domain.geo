// Gmsh project
// problem_domain_khabakhpasheva_2D.geo
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// --------------------
// Parameters
// --------------------
Lb = 12.5;          // beam length [m]
Ld = Lb;            // damping zone length [m]
Lf = 25.0;          // central fluid region length [m]
L  = 2*Ld + Lf;     // total tank length [m]
h  = 1.1;           // water depth [m]
lc = 0.1;           // target mesh size 

// Beam coordinates (centered in the middle region)
xL = Ld + (Lf - Lb)/2;   // 18.75
xR = xL + Lb;            // 31.25

// Coordinates in (x,z) with free surface at z=0 and bottom at z=-h
Point(1) = {0,  -h, 0, lc};
Point(2) = {L,  -h, 0, lc};
Point(3) = {L,   0, 0, lc};
Point(4) = {0,   0, 0, lc};
Point(5) = {xL,  0, 0, lc};
Point(6) = {xR,  0, 0, lc};

// --------------------
// Boundary curves (counterclockwise)
// --------------------
Line(1) = {1,2}; // Bottom
Line(2) = {2,3}; // Outlet

Line(3) = {3,6}; // FreeSurface right
Line(5) = {6,5}; // Structure
Line(6) = {5,4}; // FreeSurface left

Line(4) = {4,1}; // Inlet

// --------------------
// Surface
// --------------------
Curve Loop(10) = {1,2,3,5,6,4};
Plane Surface(11) = {10};

// --------------------
// Physical groups (TAGS)
// --------------------
Physical Curve("Bottom")      = {1};
Physical Curve("Outlet")      = {2};
Physical Curve("FreeSurface") = {3,6};
Physical Curve("Inlet")       = {4};
Physical Curve("Structure")   = {5};

Physical Surface("Fluid")     = {11};