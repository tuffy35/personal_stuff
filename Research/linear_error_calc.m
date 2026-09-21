clc, clear
%% mesh 3

XY = readmatrix('mesh3_coords.txt');
U = readmatrix('mesh3_linear_solution.txt');

x = XY(:,1);
y = XY(:,2);

max_error3 = 0;
for i = 1:length(XY)
    exact_sol = x(i)+y(i);
    model_sol = U(i);
    point_error = abs(exact_sol-model_sol);
    if point_error > max_error3
        max_error3 = point_error;
    end
end