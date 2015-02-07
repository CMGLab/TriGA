function [NODE,IEN] = xmesh2d(filename,varargin)
%-------------------------------------------------------------------------%
% XMESH2D stands for EXACT MESHER 2d. Given an
% arbitrary collection of NURBS curves (up through polynomial degree  3)
% in 2D, XMESH2d will generate a
% higher order triangular mesh that exactly matches the inputted geometry.

% Input:



% Output:


%-------------------------------------------------------------------------%

if nargin == 1
    options.output = false;
    options.debug = false;
elseif nargin == 2;
    options = varargin{1};
end

addpath('./Mesh2D')
% Input Parser
[P,KV,bflag,bc] = splineFileIn(filename);

% Subdividing the inputted NURBS curves into polygons
thresh = 1.01;
[node,edge,kvloc] = NURBS2poly(P,KV,thresh);

% Feeding the polygons into mesh2d
[pts,tri,kvloc] = mesh2dcurved(node,edge,P,KV,kvloc,[],options);

% Go back and perform knot insertion along the boundary at the locations
% specified by kvloc.
[bNode, bNodeflag] = meshBoundary(P,KV,bflag,kvloc);

% Go back and create higher order triangles based off of the linears
% created by mesh2d.
[node, BFLAG] = elevateMesh(pts,tri,bNode,bNodeflag,bc);

% Generate the global NODE and connectivity arrays.
[NODE,IEN] = gen_arrays(node);

% Displaying the mesh.
if options.output
    showMesh(NODE,IEN);
end

%Write out the mesh to a gambit .neu file.
gambitFileOut(filename,NODE,IEN,BFLAG)
return

function [node,BFLAG]= elevateMesh(pts,tri,bNode,bflag,bc)

side = [1 2; 2 3; 3 1];
side10 = [1 4 5 2; 2 6 7 3; 3 8 9 1];
node = cell(length(tri),1);
BFLAG = zeros(numel(bNode),4);
ctr = 1;
for ee = 1:length(tri)
    vert = pts(tri(ee,:),:);
    node{ee} = gen_net(vert);
    
    %Check to see if the current triangle is a boundary triangle.
    for bb = 1:numel(bNode)
        for ss = 1:3
            if all(single(node{ee}(side(ss,1),1:2)) == single(bNode{bb}(1,1:2))) && ...
                    all(single(node{ee}(side(ss,2),1:2)) == single(bNode{bb}(4,1:2)));
                node{ee}(side10(ss,:),:) = bNode{bb};
                
                BFLAG(ctr,1) = ee;
                BFLAG(ctr,2) = ss;
                BFLAG(ctr,3) = bflag(bb);
                BFLAG(ctr,4) = bc(bflag(bb),2);
                ctr=ctr+1;
            elseif all(single(node{ee}(side(ss,2),1:2)) == single(bNode{bb}(1,1:2))) && ...
                    all(single(node{ee}(side(ss,1),1:2)) == single(bNode{bb}(4,1:2)));
                node{ee}(side10(ss,:),:) = flipud(bNode{bb});
                BFLAG(ctr,1) = ee;
                BFLAG(ctr,2) = ss;
                BFLAG(ctr,3) = bflag(bb);
                BFLAG(ctr,4) = bc(bflag(bb),2);
                ctr=ctr+1;
            end
        end
        
    end
end

return

function [P,KV,BFLAG,BC] = splineFileIn(filename)

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


for cc = 1:nCurves
    line = fgetl(fileID);
    dims = fscanf(fileID,'%u');
    nCP = dims(1);
    lKV = dims(2);
    p   = dims(3);
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
    
for bb = 1:NBC
    line = fgetl(fileID);
    BC(bb,:) = sscanf(line,'%u');
    
end

fclose(fileID);
return
