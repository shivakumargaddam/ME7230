// DEFINED VARIABLES
hmax = .05e-3;           // mesh size -- l*0.2
crackwidth = 0.001e-3;


// POINTS
Point(1) = {0, 0, 0, hmax};
Point(2) = {1e-3, 0, 0, hmax};
Point(3) = {1e-3, 1e-3, 0, hmax};
Point(4) = {0, 1e-3, 0, hmax};
Point(5) = {0,5e-4 + crackwidth,0,hmax};
Point(6) = {5e-4,5e-4,0,hmax};
Point(7) = {0,5e-4 - crackwidth,0,hmax};
// Point(6) = {5e-4,5e-4 + crackwidth,0,hmax};
// Point(7) = {5e-4,5e-4 - crackwidth,0,hmax};
// Point(8) = {0,5e-4 - crackwidth,0,hmax};
// Point(7) = {1e-3,5e-4,hmax};


// CURVES
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
// Line(8) = {8, 1};


// CURVE LOOPS & SURFACES
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};  

Field[1] = Distance;
Field[1].CurvesList = {2};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hmax / 10;
Field[2].SizeMax = hmax;
Field[2].DistMin = 0.52e-3;
Field[2].DistMax = 1e-3;

Field[6] = Box;
Field[6].VIn = hmax / 10;
Field[6].VOut = hmax;
Field[6].XMin = 0.47e-3;
Field[6].XMax = 1e-3;
Field[6].YMin = 0e-3;
Field[6].YMax = 0.58e-3;
Field[6].Thickness = 0.3e-3;

Background Field = 6;


// PHYSICAL GROUPS
Physical Surface("domain", 1) = {1};
Physical Curve("fixed", 1) = {1};
Physical Curve("moving", 2) = {3};
// Physical Curve("crack", 4) = {5};

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;


// CREATE A MESH AND SAVE IT TO .m
Mesh 2;
// RecombineMesh;
Mesh.RecombineAll = 0;
Save "shearbox_l01_tri.m";
