%% Matlab mesh
%% square_for_reals, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 13;
msh.POS = [
0 0 0;
1 0 0;
0 1 0;
1 1 0;
0 0.5000000000013305 0;
0.4999999999986718 0 0;
1 0.4999999999986718 0;
0.5000000000013305 1 0;
0.5 0.5 0;
0.75 0.25 0;
0.25 0.25 0;
0.25 0.75 0;
0.75 0.75 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 3 5 0
 5 1 0
 1 6 0
 6 2 0
 2 7 0
 7 4 0
 4 8 0
 8 3 0
];
msh.TRIANGLES =[
 1 6 11 0
 6 10 11 0
 6 2 10 0
 11 10 9 0
 3 5 12 0
 5 11 12 0
 5 1 11 0
 12 11 9 0
 2 7 10 0
 7 13 10 0
 7 4 13 0
 10 13 9 0
 4 8 13 0
 8 12 13 0
 8 3 12 0
 13 12 9 0
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

writematrix(connectivity,'mesh1_connectivity.txt','Delimiter',' ')
writematrix(coords,'mesh1_coords.txt','Delimiter',' ')
