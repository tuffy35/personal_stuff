clc,clear,close all
% mesh 0
nodes0 = load('mesh0_coords.txt');        % n x 2 array: [x y]
triangles0 = load('mesh0_connectivity.txt'); % m x 3 array: indices into nodes

% Separate x and y coordinates
x0 = nodes0(:,1) ;
y0 = nodes0(:,2) ;

triangles0 = triangles0 + 1;

% mesh 2
nodes2 = load('mesh2_coords.txt');        % n x 2 array: [x y]
triangles2 = load('mesh2_connectivity.txt'); % m x 3 array: indices into nodes

% Separate x and y coordinates
x2 = nodes2(:,1) ;
y2 = nodes2(:,2) ;

triangles2 = triangles2 + 1;

% mesh 4

nodes4 = load('mesh4_coords.txt');        % n x 2 array: [x y]
triangles4 = load('mesh4_connectivity.txt'); % m x 3 array: indices into nodes

% Separate x and y coordinates
x4 = nodes4(:,1) ;
y4 = nodes4(:,2) ;

triangles4 = triangles4 + 1;

% mesh 5
nodes6 = load('mesh5_coords.txt');        % n x 2 array: [x y]
triangles5 = load('mesh5_connectivity.txt'); % m x 3 array: indices into nodes

% Separate x and y coordinates
x5 = nodes6(:,1) ;
y5 = nodes6(:,2) ;

triangles5 = triangles5 + 1;

%trimeshImproved(triangles0,x0,y0,1)

%trimeshImproved(triangles2,x2,y2,2)

%trimeshImproved(triangles4,x4,y4,3)

trimeshImproved(triangles5,x5,y5,1)

