clc
clear
close all

% ------------------------------------------------------------------------%
% This script uses the method of manufactured solutions to validate the
% convergence rate of numerical methods.

% The problem geometry is a 8x8 square centered at the origin, with a unit
% circle hole cut in the center. 

% The manufactured solution is formed so that the temperature is zero on
% all of the boundaries. The sides of the square are set as diriclet
% boundary conditions and the edge of the circle is set as a Neumann
% boundary condition. 
% ------------------------------------------------------------------------%

addpath('xmesh')
addpath('xmesh/Mesh2D')
filename = 'plate';
xmesh2d(filename);

% Define the function for the manufactured solution.
phi = @(x,y) -(1-sqrt(x^2+y^2))*(x+4)*(x-4)*(y+4)*(y-4);

gradPhi = @(x,y) [(x*(-4+y)*(4+y)*(3*x^2-2*(8-y^2+sqrt(x^2+y^2))))/sqrt(x^2+y^2),...
    ((-4+x)*(4+x)*y*(2*x^2+3*y^2-2*(8+sqrt(x^2+y^2))))/sqrt(x^2+y^2)];
heatGen = @(x,y)  -(2 *x^4+x^2 *(13 *y^2-2 *(sqrt(x^2+y^2)+72))+2 *(-y^2 *(sqrt(x^2+y^2)+72)+32 *(sqrt(x^2+y^2)+4)+y^4))/sqrt(x^2+y^2);



[NODE,IEN,temp] = heat2d(filename,heatGen,[0 0]);

showResults(filename)

% Method of manufactured solutions
[L2, H1] = calcNorm(filename,phi,gradPhi);


nref = 5;
q = zeros(nref+1,2);
qtemp = quality(NODE(:,1:2),IEN(1:3,:)');
q(1,:) = [mean(qtemp),min(qtemp)];
% Refinement
for i = 1:nref
    refineMesh(filename);
    filename = [filename,'ref']; %#ok<AGROW>
    [NODE,IEN,temp] = heat2d(filename,heatGen,[0 0 0]);
    qtemp = quality(NODE(:,1:2),IEN(1:3,:)');
    q(i+1,:) = [mean(qtemp),min(qtemp)];
    [L2(i+1), H1(i+1)] = calcNorm(filename,phi,gradPhi);
end