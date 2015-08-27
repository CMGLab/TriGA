function [NODE,IEN] = xmesh2d(filename,varargin)
%-------------------------------------------------------------------------%
% XMESH2D stands for EXACT MESHER 2d. Given an arbitrary collection of
% NURBS curves up through polynomial degree  3 in 2D, XMESH2d will generate
% a higher order triangular mesh that exactly matches the input geometry.
%
% INPUT:
% filename: The filename of a .spline file containing the curves that
% define the problem geometry.
%
% Output:
% NODE: A nNodesx3 array containing the coordinates and weights of the
% nodes that make up the mesh.
%
% IEN: a 10xnel array containing the mesh connectivity information.
%
% filename.neu: xmesh2d also writes out the mesh information int he
% standard gambit neu trail file format to the file <filename>.neu
%-------------------------------------------------------------------------%

if nargin == 1
    options.output = false;
    options.debug = false;
    options.smooth = true;
elseif nargin == 2;
    options = varargin{1};
    if ~exist('options.output') %#ok<*EXIST>
        options.output = false;
    end
    if ~exist('options.debug')
        options.debug = false;
    end
    if ~exist('options.smooth')
        options.smooth = true;
    end
    
    
end

addpath('~/TriGA/xmesh/Mesh2D')
addpath('./xmesh/Mesh2D')
addpath('./Mesh2D')
addpath('.\xmesh\Mesh2D')
addpath('.\Mesh2D')
addpath('..\xmesh\Mesh2D')
addpath('..\Mesh2D')
clc

% Input Parser
[P,KV,p,bflag,bc,face] = splineFileIn(filename);

% Normalize the knot vectors to span [0 1]
for kk = 1:numel(KV)
    KV{kk} = KV{kk}/KV{kk}(end);
end

% Subdividing the inputted NURBS curves into polygons
thresh = 1.01;
[node,edge,kvloc] = NURBS2poly(P,KV,thresh);

% Feeding the polygons into mesh2d
[pts,tri,kvloc] = xmeshfaces(node,edge,P,KV,kvloc,face,[],options);

% Go back and perform knot insertion along the boundary at the locations
% specified by kvloc.
[bNode, bNodeflag] = meshBoundary(P,KV,bflag,kvloc);

% Go back and create higher order triangles based off of the linears
% created by mesh2d.
[node,linnode,BFLAG,CFLAG] = elevateMesh(pts,tri,bNode,bNodeflag,bc,p);

% Generate the global NODE and connectivity arrays.
[NODE,IEN] = gen_arrays(node);
[linNODE,~] = gen_arrays(linnode);

% Displaying the mesh.
if options.output
    showMesh(NODE,IEN);
end

if options.smooth
    % Solve a Laplace problem to smooth the weights.
    smoothNODE = smoothWeights(NODE, IEN, BFLAG);
    smoothNODE = smoothNODE(:,3);
end
NODE(:,3) = smoothNODE;
%Write out the mesh to a gambit .neu file.
gambitFileOut(filename,NODE,IEN,BFLAG,CFLAG)
return

function [node,linnode,BFLAG,CFLAG]= elevateMesh(pts,tri,bNode,bflag,bc,p)

side = [1 2; 2 3; 3 1];
side10 = [1 4 5 2; 2 6 7 3; 3 8 9 1];
node = cell(length(tri),1);
linnode = cell(length(tri),1);

BFLAG = zeros(numel(bNode),4);
ctr = 1;
d = 14;

% Define CFLAG as a flag if the current element is curved.
CFLAG = false(length(tri),1);
for ee = 1:length(tri)
    vert = pts(tri(ee,:),:);
    linnode{ee} = round(gen_net(vert)*10^d)/10^d;
    node{ee} = round(gen_net(vert)*10^d)/10^d;
    %Check to see if the current triangle is a boundary triangle.
    for bb = 1:numel(bNode)
        
        % Check to see if the triangle has a side on the boundary.
        for ss = 1:3
            if all(single(node{ee}(side(ss,1),1:2)) == single(bNode{bb}(1,1:2))) && ...
                    all(single(node{ee}(side(ss,2),1:2)) == single(bNode{bb}(4,1:2)));
                node{ee}(side10(ss,:),:) = bNode{bb};
                linnode{ee}(side10(ss,:),3) = bNode{bb}(:,3);
                
                BFLAG(ctr,1) = ee;
                BFLAG(ctr,2) = ss;
                BFLAG(ctr,3) = bflag(bb,1);
                BFLAG(ctr,4) = bc(bflag(bb,1),2);
                
                if p(bflag(bb,2))>1
                    CFLAG(ee) = true;
                end
                ctr=ctr+1;
            elseif all(single(node{ee}(side(ss,2),1:2)) == single(bNode{bb}(1,1:2))) && ...
                    all(single(node{ee}(side(ss,1),1:2)) == single(bNode{bb}(4,1:2)));
                node{ee}(side10(ss,:),:) = flipud(bNode{bb});
                linnode{ee}(side10(ss,:),3) = flipud(bNode{bb}(:,3));
                
                BFLAG(ctr,1) = ee;
                BFLAG(ctr,2) = ss;
                BFLAG(ctr,3) = bflag(bb,1);
                BFLAG(ctr,4) = bc(bflag(bb,1),2);
                
                if p(bflag(bb,2))>1
                    CFLAG(ee) = true;
                end
                ctr=ctr+1;
            end
        end
        
        % Also check to see if the triangle has a vertex on the boundary.
        % If it does, it is part of the blending layer, so it needs to be
        % flagged.
        for vv = 1:3
            if single(node{ee}(vv,1:2)) == single(bNode{bb}(1,1:2)) || ...
                    single(node{ee}(vv,1:2)) == single(bNode{bb}(4,1:2));
                if p(bflag(bb,2))>1
                    CFLAG(ee) = true;
                end
                ctr=ctr+1;
            end
        end
        
    end
end

return

function [node] = gen_net(vert)
% ----------------------------------gen_net------------------------------------%
% Generates the control net for the current triangular Bezier patch. Reads in
% the verticies of a triangle, and linearly interpolates the triangle to make a
% 10 control point bezier element.

% INPUT:
% vert: a matrix in the form [x1 y1; x2 y2; x3 y3] containing the verticies of a
%  linear triangle, ordered counter clockwise from the bottom left vertex.

%     v3
%     | \
%     |  \
%     |   \
%     |    \
%     v1----v2

% OUTPUT:
% node: a 10x3 array containing the control point coordinates and weights.
% The standard node ordering is shown below.
%
%     3
%     |\
%     | \
%     8  7
%     |   \
%     |    \
%     9 10  6
%     |      \
%     |       \
%     1--4--5--2
%
%------------------------------------------------------------------------------%


node(:,1:2) = [ vert(1,1:2);...
    vert(2,1:2);...
    vert(3,1:2);...
    vert(1,1:2) + (vert(2,1:2)-vert(1,1:2))/3;...
    vert(1,1:2) + (vert(2,1:2)-vert(1,1:2))*2/3;...
    vert(2,1:2) + (vert(3,1:2)-vert(2,1:2))/3;...
    vert(2,1:2) + (vert(3,1:2)-vert(2,1:2))*2/3;...
    vert(3,1:2) + (vert(1,1:2)-vert(3,1:2))/3;...
    vert(3,1:2) + (vert(1,1:2)-vert(3,1:2))*2/3];

node(10,1:2) = node(9,1:2) + (node(6,1:2)-node(9,1:2))/2;
if size(vert,2) ==3
    node(4:10,3) = ones(7,1);
    node(1:3,3) = vert(1:3,3);
else
    node(:,3) = ones(10,1);
end

return
function [P,KV,p,BFLAG,BC,FACE] = splineFileIn(filename)

filename = [filename,'.spline'];
fileID = fopen(filename,'r');

% Read in the file header
for i=1:2
    line = fgetl(fileID);
end

% Find the total number of curves
line = fgetl(fileID);
dims = sscanf(line, '%f');
nCurves = dims(1);

KV = cell(1,nCurves);
P  = cell(1,nCurves);
BFLAG = cell(1,nCurves);
p = zeros(1,nCurves);

for cc = 1:nCurves
    line = fgetl(fileID); %#ok<*NASGU>
    dims = fscanf(fileID,'%u');
    nCP = dims(1);
    lKV = dims(2);
    p(cc)   = dims(3);
    line=fgetl(fileID);
    
    for ii = 1:nCP
        line = fgetl(fileID);
        tmpx = sscanf(line, '%lf');
        P{cc}(ii,1) = tmpx(1); P{cc}(ii,2) = tmpx(2); P{cc}(ii,3) = tmpx(3);
    end
    line = fgetl(fileID);
    line = fgetl(fileID);
    KV{cc} = sscanf(line, '%lf');
    KV{cc} = KV{cc}';
    line = fgetl(fileID);
    nKspan = length(unique(KV{cc}))-1;
    
    for bb = 1:nKspan
        line = fgetl(fileID);
        BFLAG{cc}(bb,:) = sscanf(line, '%u');
    end
    
end
line = fgetl(fileID);line = fgetl(fileID);line = fgetl(fileID);
NBC = sscanf(line, '%u');

BC = zeros(NBC,2);
for bb = 1:NBC
    line = fgetl(fileID);
    BC(bb,:) = sscanf(line,'%u');
    
end

line = fgetl(fileID);

if line ~= -1
    fgetl(fileID);
    
    line = fgetl(fileID);
    nFace = sscanf(line,'%u');
    FACE = cell(1,nFace);
    for ff = 1:nFace
        fgetl(fileID);
        line = fgetl(fileID);
        nCurve = sscanf(line,'%u');
        for cc = 1:nCurve
            line = fgetl(fileID);
            FACE{ff}(cc) = sscanf(line,'%u');
        end
    end
else
    FACE{1} = 1:nCurves;
end

fclose(fileID);
return
