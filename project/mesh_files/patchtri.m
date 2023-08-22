%  Matlab mesh
% R, Created by Gmsh
% ASCII
clear msh;
msh.nbNod = 4;
msh.POS = [
0 0 0;
1 0 0;
1 1 0;
0 1 0
]*1e-3;
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES = [
1 2 1;
3 4 2;
1 2 1;
3 4 2;
1 2 1;
3 4 2;
1 2 1;
3 4 2
];
msh.TRIANGLES =[
1 2 3 0;
1 3 4 0
];

