% ------------------------------------------------------------------------%
% TRIGADEMO
% This script demonstrates the use of TriGA over two different geometries.
% It also uses the method of manufatured solutions to verify that the
% heat2d solver does indeed converge as expected. 
%
% NOTICE! This can take a while to run! 
%
% TriGA has not been optimized for speed yet, and is very innefficient.
% Takes about 5 minutes to run on a quad core MacBook Pro. 
% ------------------------------------------------------------------------%
tic
clc
clear
close all
addpath('xmesh')
addpath('xmesh/Mesh2D')

% ------------------------------------------------------------------------%
% EXAMPLE 1: PLATE AND HOLE
% This problem consists of a 8x8 plate centered at the origin with a unit
% circle hole on the middle. Internal heat generation is defined based off
% of a manufactured solution, and all boundary conditions have homogeneous
% diriclet boundary conditions.
% ------------------------------------------------------------------------%

filename = './plateandhole/plateandhole';

% Define the function for the manufactured solution.
phi = @(x,y) -(1-sqrt(x^2+y^2))*(x+4)*(x-4)*(y+4)*(y-4);

gradPhi = @(x,y)...
    [(x*(-4+y)*(4+y)*(3*x^2-2*(8-y^2+sqrt(x^2+y^2))))/sqrt(x^2+y^2),...
    ((-4+x)*(4+x)*y*(2*x^2+3*y^2-2*(8+sqrt(x^2+y^2))))/sqrt(x^2+y^2)];
heatGen = @(x,y) -(2 *x^4+x^2 *(13 *y^2-2 *(sqrt(x^2+y^2)+72))+2 *...
    (-y^2 *(sqrt(x^2+y^2)+72)+32 *(sqrt(x^2+y^2)+4)+y^4))/sqrt(x^2+y^2);

BC = [0 0 0];

D = eye(2);

% Generating the mesh.
xmesh2d(filename);

% Calling heat2d
heat2d(filename,heatGen,BC,D);

% Visualizing results
showResults(filename,1/10)

% Calculating the residual.
nref = 1;
[L2, H1] = calcNorm(filename,phi,gradPhi);
h = zeros(1,nref+1);
h(1) = 1;

% Refinement
for i = 1:nref
    h(i+1) = h(i)/2;
    refineMesh(filename);
    filename = [filename,'ref']; %#ok<AGROW>
    heat2d(filename,heatGen,BC,D);
    [L2(i+1), H1(i+1)] = calcNorm(filename,phi,gradPhi);
end

% Plotting Convergence
figure
loglog(h,L2,'or')
hold on
loglog(h,L2,'--r')
set(gca,'xdir','reverse','FontSize',16)
xlabel('Mesh Size')
ylabel('L2 Norm')
title('Convergence for Example 1')

figure
loglog(h,H1,'or')
hold on
loglog(h,H1,'r--')
set(gca,'xdir','reverse','FontSize',16)
xlabel('Mesh Size')
ylabel('H1 Norm')
title('Convergence for Example 1')

%% ------------------------------------------------------------------------%
% EXAMPLE 2: QUARTER ANNULUS
% This problem consists of a quarter annulus with r_i = 1 and r_o = 2. A
% constant Neumann boundary condition is applied at the inner radius and a
% homogeneous diriclet boundary condition is applied at the outer radius.
% No-flux boundary conditions are applied at the ends of the annulus. 
% ------------------------------------------------------------------------%
clear
filename = './quarterannulus/quarterannulus';

r = @(x,y) sqrt(x^2+y^2);
phi = @(x,y) -log(r(x,y)/2)/log(2);
gradPhi = @(x,y) [-2*x/(log(4)*(x^2+y^2)),-2*y/(log(4)*(x^2+y^2))];
heatGen  = @(x,y) 0;
BC = [0 0 2*0.721347520444482];
D = eye(2);

% Generating the mesh.
xmesh2d(filename);

% Calling heat2d
[NODE,IEN,temp] = heat2d(filename,heatGen,BC,D);

% Visualizing results
showResults(filename,1/30)

nref = 3;
L2 = zeros(1,nref+1); %#ok<*PREALL>
H1 = zeros(1,nref+1);
h = zeros(1,nref+1);
h(1) = 1;

% Calculating the residual.
[L2, H1] = calcNorm(filename,phi,gradPhi);
% Refinement
for i = 1:nref
    h(i+1) = h(i)/2;
    refineMesh(filename);
    filename = [filename,'ref']; %#ok<AGROW>
    heat2d(filename,heatGen,BC,D);
    [L2(i+1), H1(i+1)] = calcNorm(filename,phi,gradPhi);
end

% Plotting Convergence
figure
loglog(h,L2,'or')
hold on
loglog(h,L2,'--r')
set(gca,'xdir','reverse','FontSize',16)
xlabel('Mesh Size')
ylabel('L2 Norm')
title('Convergence for Example 2')

figure
loglog(h,H1,'or')
hold on
loglog(h,H1,'r--')
set(gca,'xdir','reverse','FontSize',16)
xlabel('Mesh Size')
ylabel('H1 Norm')
title('Convergence for Example 2')

toc
