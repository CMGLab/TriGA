clc
clear
close all

filename = './circle/circA01curved';
% 
% phi = @(x,y) 1-exp(-(sqrt(x^2+y^2)-1)^2);
% gradPhi = @(x,y) [2*exp(-(-1+sqrt(x^2+y^2))^2)*x*(-1+sqrt(x^2+y^2)),...
%                     2*exp(-(-1+sqrt(x^2+y^2))^2)*y*(-1+sqrt(x^2+y^2))]/sqrt(x^2+y^2);
% heatGen = @(x,y) exp(-(sqrt(x^2+y^2)-1)^2)*(-4*x^2*(sqrt(x^2+y^2)-2)-4*y^2*(sqrt(x^2+y^2)-2)-2)/sqrt(x^2+y^2);
phi = @(x,y) -(1-sqrt(x^2+y^2));
gradPhi = @(x,y) [x,y]/sqrt(x^2+y^2);
heatGen = @(x,y) -1/sqrt(x^2+y^2);


BC = 0;
D = eye(2);

heat2d(filename,heatGen,BC,D);

% Calculating the residual.
nref = 3;
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