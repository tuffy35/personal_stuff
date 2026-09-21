%% EXPONENTIAL
exp_error = [0.5	0.003097
0.25	0.000792
0.125	2.14E-04
0.0625	5.00E-05
0.03125	1.28E-05
0.015625	3.30E-06
0.0078125	8.59E-07];

figure(1)
loglog(exp_error(:,1),exp_error(:,2))
%title('2D Mesh Refinement Study: Exponential')
ylabel('max error')
xlabel('mesh spacing')
grid on

%% QUADRATIC
quad_error = [0.5	0.549999999
0.25	0.140624999
0.125	0.037764442
0.0625	0.008997545
0.03125	0.002300076
0.015625	0.000595328
0.0078125	0.000155095];

figure(2)
loglog(quad_error(:,1),quad_error(:,2))
%title('2D Mesh Refinement Study: Quadratic')
ylabel('max error')
xlabel('mesh spacing')
grid on

%% Sinusoid
sin_error = [0.25	0.999999999
0.125	0.57253065
0.0625	0.18333084
0.03125	0.03544737
0.015625	0.00868873
0.0078125	0.00285811];

figure(3)
loglog(sin_error(:,1),sin_error(:,2))
%title('2D Mesh Refinement Study: Sinusoidal')
ylabel('max error')
xlabel('mesh spacing')
grid on

%% Cubic
cubic_error = [0.5	0.003097132
0.25	0.000792252
0.125	0.000214376
0.0625	5.00E-05
0.03125	1.28E-05
0.015625	3.30E-06
0.0078125	8.25E-07];

figure(4)
loglog(cubic_error(:,1),cubic_error(:,2))
%title('2D Mesh Refinement Study: Cubic')
ylabel('max error')
xlabel('mesh spacing')
grid on

%% Combined

subplot(2,2,1)
loglog(exp_error(:,1),exp_error(:,2),'-s')
title('2D Mesh Refinement Study: Exponential')
ylabel('max error')
xlabel('mesh spacing')
grid on

subplot(2,2,2)
loglog(quad_error(:,1),quad_error(:,2),'-s')
title('2D Mesh Refinement Study: Quadratic')
ylabel('max error')
xlabel('mesh spacing')
grid on

subplot(2,2,3)
loglog(sin_error(:,1),sin_error(:,2),'-s')
title('2D Mesh Refinement Study: Sinusoidal')
ylabel('max error')
xlabel('mesh spacing')
grid on

subplot(2,2,4)
loglog(cubic_error(:,1),cubic_error(:,2),'-s')
title('2D Mesh Refinement Study: Cubic')
ylabel('max error')
xlabel('mesh spacing')
grid on