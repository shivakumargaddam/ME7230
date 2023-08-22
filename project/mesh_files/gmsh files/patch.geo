// DEFINED VARIABLES
delta = 1.4e-5;           // mesh size -- l*0.2


// POINTS
Point(1) = {0, 0, 0, delta};
Point(2) = {1e-3, 0, 0, delta};
Point(3) = {1e-3, 1e-3, 0, delta};
Point(4) = {0, 1e-3, 0, delta};


// CURVES
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


// CURVE LOOPS & SURFACES
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};


// PHYSICAL GROUPS
Physical Surface("domain", 1) = {1};
Physical Curve("fixed", 1) = {1};
Physical Curve("moving", 2) = {3};

Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 1;

// CREATE A MESH AND SAVE IT TO .m
Mesh 2;
// RecombineMesh;
Mesh.RecombineAll = 1;
Save "patchquad3.m";
