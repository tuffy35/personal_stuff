clc, clear, close all

Sol = readmatrix('solution_exponential_8.txt');
% x^2 + 2y^2 + 3z^2 + 4w^2
XYZW = readmatrix('coordinates_mesh_8.txt');
pents = readmatrix('cells_mesh_8.txt')+1;

% x = XYZW(:,1);
% y = XYZW(:,2);
% z = XYZW(:,3);
% w = XYZW(:,4);

dx = 0.01;
dy = dx;
dz = dx;
dw = dx;

x = 0:dx:1;
y = 0:dy:1;
z = 0.5;
w = 0.5;

[X,Y] = meshgrid(x,y);

[nrows, ncols] = size(X);
U = zeros(nrows,ncols);

ethresh = 0.00000001;

for i = 1:nrows
    for j = 1:ncols
        
        xlocal = X(i,j);
        ylocal = Y(i,j);

        zlocal = z;
        wlocal = w;

        fvalue = 0;

        for k = 1:length(pents)

            ia = pents(k,1);
            ib = pents(k,2);
            ic = pents(k,3);
            id = pents(k,4);
            ie = pents(k,5);

            va = XYZW(ia,:);
            vb = XYZW(ib,:);
            vc = XYZW(ic,:);
            vd = XYZW(id,:);
            ve = XYZW(ie,:);

            column_va = [va 1]';
            column_vb = [vb 1]';
            column_vc = [vc 1]';
            column_vd = [vd 1]';
            column_ve = [ve 1]';

            Mat = [column_va, column_vb, column_vc, column_vd, column_ve];
            b = [xlocal ylocal zlocal wlocal 1]';

            bary = Mat\b;

            bary_flag = 1;
            for m = 1:length(bary)

                if((bary(m) < - ethresh) || (bary(m) > (1 + ethresh)) )
                    bary_flag = 0;
                end
            end

            if(bary_flag == 1)

                fvalue = bary(1)*Sol(ia) + bary(2)*Sol(ib) + bary(3)*Sol(ic) + bary(4)*Sol(id) + bary(5)*Sol(ie);
                break;
            end
        end

        U(i,j) = fvalue;
    end
end
%% Plotting

% camera angle
a1 = -14.0292;
a2 = 6.1788;

figure(1)
surf(X,Y,U,'EdgeColor','none')
axis square
colorbar
colormap(turbo)

% 
% Uexact = X.^2 + 2*Y.^2 + 3*z^2 + 4*w^2;
% 
% figure(2)
% surf(X,Y,Uexact,'EdgeColor','none')
% axis square
% colorbar
% colormap(turbo)