SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 2.2;

// --------------------------------------------------
// Pipe diameter
// --------------------------------------------------
D  = 0.40;                 // diameter

// --------------------------------------------------
// Domain parameters
// --------------------------------------------------
Ld = 40.0;      // damping zone length [m]
Lf = 35.0;      // central fluid region length [m]
L  = 2*Ld + Lf; // total tank length [m] = 75.0
h  = 19.0;      // water depth [m]
lc = D / 5 ;       // mesh size near pipe

// Damping zone boundaries
xdin  = Ld;         // 20.0
xdout = Ld + Lf;    // 55.0

// --------------------------------------------------
// Pipe parameters
// --------------------------------------------------
R  = D/2;                 // radius
sub = 2.0*D/3.0;          // submerged depth
zc = -(sub - R);          // pipe center z-coordinate
xc = L/2;                 // pipe center x-coordinate

xoff = Sqrt(R*R - zc*zc); // half waterline chord
xL = xc - xoff;           // left waterline intersection
xR = xc + xoff;           // right waterline intersection

// --------------------------------------------------
// Outer tank boundary
// --------------------------------------------------
Point(1) = {0,    -h, 0, 0.8};
Point(2) = {L,    -h, 0, 0.8};
Point(3) = {L,     0, 0, 0.8};
Point(4) = {0,     0, 0, 0.8};

// Damping zone / free surface transitions
Point(9)  = {xdin,  0, 0, lc};
Point(10) = {xdout, 0, 0, lc};

// free-surface intersection points
Point(5) = {xL, 0, 0, lc};
Point(6) = {xR, 0, 0, lc};

// circle center
Point(7) = {xc, zc, 0, lc};

// lowest point of circle
Point(8) = {xc, zc - R, 0, lc};

// --------------------------------------------------
// Boundary lines
// --------------------------------------------------
Line(1)  = {1,2};    // Bottom
Line(2)  = {2,3};    // Outlet
Line(3)  = {4,1};    // Inlet

Line(11) = {3,10};   // FreeSurface outlet damping
Line(4)  = {10,6};   // FreeSurface right of pipe
Line(5)  = {5,9};    // FreeSurface left of pipe
Line(12) = {9,4};    // FreeSurface inlet damping

// Wetted pipe boundary: lower semicircle-like arc from right to left
Circle(6) = {6,7,8};
Circle(7) = {8,7,5};

// --------------------------------------------------
// Fluid surface
// --------------------------------------------------
Curve Loop(10) = {1,2,11,-4,-6,-7,-5,12,3};
Plane Surface(11) = {10};

// --------------------------------------------------
// Physical groups
// --------------------------------------------------
Physical Curve("Bottom")      = {1};
Physical Curve("Outlet")      = {2};
Physical Curve("Inlet")       = {3};
Physical Curve("FreeSurface") = {11,4,5,12};
Physical Curve("Pipe")        = {6,7};

Physical Surface("Fluid")     = {11};

Physical Point("BottomLeft")    = {1};
Physical Point("DampingInlet")  = {9};
Physical Point("DampingOutlet") = {10};