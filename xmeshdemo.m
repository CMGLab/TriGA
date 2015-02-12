% ------------------------------------------------------------------------%
% XMESHDEMO This script demonsates the capabilites of xmesh. 
% ------------------------------------------------------------------------%
clc
clear
close all
addpath xmesh
addpath xmesh/Mesh2D

% ------------------------------------------------------------------------%
% EXAMPLE 1 is just the unit circle.
% ------------------------------------------------------------------------%
in = input('Example 1 is the unit circle, would you like to see it? (y/n): ','s');
if strcmp(in,'y')
    filename = './circle/circle';
    xmesh2d(filename);
    showMesh(filename)
    
    in = input('xmesh can also do refinement. Would you like to see an example? (y/n): ','s');
    if strcmp(in,'y')
        refineMesh(filename);
        showMesh([filename,'ref'])
    end
end

% ------------------------------------------------------------------------%
% EXAMPLE 2 is a plate and hole.
% ------------------------------------------------------------------------%
in = input('Example 2 is a plate with a hole, would you like to see it? (y/n): ','s');
if strcmp(in,'y')
    close all
    filename = './plateandhole/plateandhole';
    xmesh2d(filename);
    showMesh(filename)
end


% ------------------------------------------------------------------------%
% EXAMPLE 3 is a nozzle with a tight constriction.
% ------------------------------------------------------------------------%
in = input('EXAMPLE 3 is a nozzle with a tight constriction, would you like to see it? (y/n): ','s');
if strcmp(in,'y')
    close all
    filename = './nozzle/nozzle';
    xmesh2d(filename);
    showMesh(filename)
end

% ------------------------------------------------------------------------%
% EXAMPLE 4 is a flange.
% ------------------------------------------------------------------------%
in = input('EXAMPLE 4 is flange, would you like to see it? (y/n): ','s');
if strcmp(in,'y')
    close all
    filename = './flange/flange';
    xmesh2d(filename);
    showMesh(filename)
end

% ------------------------------------------------------------------------%
% EXAMPLE 5 is a nozzle with a tight constriction.
% ------------------------------------------------------------------------%
in = input('EXAMPLE 5 is an airfoil, would you like to see it? (y/n): ','s');
if strcmp(in,'y')
    close all
    filename = './airfoil/airfoil';
    xmesh2d(filename);
    showMesh(filename)
end