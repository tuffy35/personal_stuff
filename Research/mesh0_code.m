%% Matlab mesh
%% square_for_reals, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 5;
msh.POS = [
0 0 0;
1 0 0;
0 1 0;
1 1 0;
0.5 0.5 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 3 1 0
 1 2 0
 2 4 0
 4 3 0
];
msh.TRIANGLES =[
 1 2 5 0
 3 1 5 0
 2 4 5 0
 4 3 5 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
];

msh.TRIANGLES = msh.TRIANGLES(:,1:3);
connectivity = msh.TRIANGLES-1;
coords = msh.POS;

writematrix(connectivity,'mesh0_connectivity.txt','Delimiter',' ')
writematrix(coords,'mesh0_coords.txt','Delimiter',' ')
