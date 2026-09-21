%% mesh 0
clc, clear
XY = readmatrix('mesh0_coords.txt');
U = readmatrix('mesh0_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error0 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error0
        max_error0 = point_error;
    end
end

%% mesh 1
XY = readmatrix('mesh1_coords.txt');
U = readmatrix('mesh1_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error1 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error1
        max_error1 = point_error;
    end
end

%% mesh 2

XY = readmatrix('mesh2_coords.txt');
U = readmatrix('mesh2_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error2 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error2
        max_error2 = point_error;
    end
end
%% mesh 3

XY = readmatrix('mesh3_coords.txt');
U = readmatrix('mesh3_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error3 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error3
        max_error3 = point_error;
    end
end
%% mesh 4

XY = readmatrix('mesh4_coords.txt');
U = readmatrix('mesh4_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error4 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error4
        max_error4 = point_error;
    end
end
%% mesh 5

XY = readmatrix('mesh5_coords.txt');
U = readmatrix('mesh5_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error5 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error5
        max_error5 = point_error;
    end
end
%% mesh 6

XY = readmatrix('mesh6_coords.txt');
U = readmatrix('mesh6_cubic_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error6 = 0;
for i = 1:length(XY)
    exact_sol = x(i)*x(i)*y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error6
        max_error6 = point_error;
    end
end
%%
% values from exponential values calculated in python
max_error0 = 2;
max_error1 = 1;
max_error2 = 0.57253065;
max_error3 = 0.18333084;
max_error4 = 0.03544737;
max_error5 = 0.00868873;
max_error6 = 0.00285811;

max_error = [max_error0, max_error1 ,max_error2,max_error3,max_error4,max_error5,max_error6];
refinement = [0,1,2,3,4,5,6];

h_size = zeros(size(max_error));
for i = 1:length(h_size)
    h_size(i) = 1/2^i;
end


figure(1)
plot(refinement,max_error)
xlabel('refinement level')
xticks([0,1,2,3,4,5,6])
ylabel('error')
title('2D Mesh Refinement Study')
grid on

figure(2)
loglog(h_size,max_error)
title('2D Mesh Refinement Study')
ylabel('max error')
xlabel('mesh spacing')
grid on

slope1 = log(max_error(2)/max_error(1))/log(2);
slope2 = log(max_error(3)/max_error(2))/log(2);
slope3 = log(max_error(4)/max_error(3))/log(2);
slope4 = log(max_error(5)/max_error(4))/log(2);
slope5 = log(max_error(6)/max_error(5))/log(2);
slope6 = log(max_error(7)/max_error(6))/log(2);
disp(slope5)

error_table = [h_size;max_error]';
%writematrix(error_table,'mesh_refinement_values_quadratic2.csv')