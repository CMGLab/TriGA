function [IEN,NODE,node,bflag] = gen_mesh(P,KV,h,options)
%------------------------------------------------------------------------------%
% Gen mesh takes as inputs n NURBS curves in the form of n knot vectors and
% n lists of control points. It also takes in h values describing the
% resolution to break up the NURBS curves with.

% Input:
% KV
% P
% h


%
% Output:
% NODE
% IEN
% bflag
%------------------------------------------------------------------------------%

addpath('~/Dropbox/TriGA2D/mesh2D/')

% Load the geometry input file

% See if the option to lift he boundary was specified.

if options.lift
    % Loop through each curve defined in the input file
    triCtr = 1;
    nodeCtr = 1;
    edgeCtr = 1;
    for nn = 1:numel(P)
        [bNodeTemp,bEdgeTemp,nodeTemp] = meshBoundary(P{nn},KV{nn},h(nn),options);
        
        for j = 1:2*h(nn)
            node{triCtr+j-1} = nodeTemp{j};
        end
        bNode(nodeCtr:nodeCtr+h(nn)-1,:) = bNodeTemp;
        bEdge(edgeCtr:edgeCtr+h(nn)-1,:) = bEdgeTemp+edgeCtr-1;
        triCtr = triCtr+2*h(nn);
        nodeCtr = nodeCtr+h(nn);
        edgeCtr = edgeCtr+h(nn);
    end
    
    bNode(nodeCtr:nodeCtr+length(PSLG)-1,:) = PSLG;
    bEdge(edgeCtr:edgeCtr+length(e)-1,:) = e+max(bEdge(:));
    
    % Call the third party mesher
    options.output  = false;
    [pts,tri] = mesh2d(bNode,bEdge,[],options);
    
    for tt = 1:length(tri)
        vert = pts(tri(tt,:),:);
        node{triCtr+tt-1} = gen_net(vert);
    end
    
else
    
    % Loop over the curves specified in the input file
    nodeCtr = 1;
    edgeCtr = 1;
    for nn = 1:numel(P)
        
        [bNodeTemp,bEdgeTemp,bPtsTemp,bflagTemp] = meshBoundary(P{nn},KV{nn},h(nn),options);
        bNode(nodeCtr:nodeCtr+h(nn)-1,:) = bNodeTemp;
        bEdge(edgeCtr:edgeCtr+h(nn)-1,:) = bEdgeTemp+edgeCtr-1;
        bPts(edgeCtr:edgeCtr+h(nn)-1) = bPtsTemp;
        bflag(edgeCtr:edgeCtr+h(nn)-1,1) = ones(size(bflagTemp))*nn;
        bflag(edgeCtr:edgeCtr+h(nn)-1,2) = bflagTemp;
        
        nodeCtr = nodeCtr+h(nn);
        edgeCtr = edgeCtr+h(nn);
        
    end
    
    % Call the third party mesher
    options.output  = false;
    [pts,tri,~] = mesh2d(bNode,bEdge,[],options);
    
    bNode (:,3) = 1:length(bNode);
    
    edgeCtr = 1;
    for nn = 1:numel(P)
        edgeTemp = zeros(h(nn),2);
        pts2bndTemp = zeros(h(nn),2);
        
        hTemp  = [0 cumsum(h)];
        bNodeTemp = bNode(hTemp(nn)+1: hTemp(nn+1),:);
        for ii = 1:length(bNodeTemp)
            jj = find(bNodeTemp(ii,1) == pts(:,1) & bNodeTemp(ii,2) == pts(:,2));
            pts2bndTemp(ii,:) = [ii,jj];
        end
        
        for ii = 1:length(pts2bndTemp)-1
            edgeTemp(ii,:) = [pts2bndTemp(ii,2), pts2bndTemp(ii+1,2)];
        end
        edgeTemp(length(pts2bndTemp),:) = [pts2bndTemp(end,2),pts2bndTemp(1,2)];
        edge(edgeCtr:edgeCtr+h(nn)-1,:) = edgeTemp;
        edgeCtr = edgeCtr+h(nn);
        
    end
    
    bflag(:,3:4) = bflag;
    bflag(:,1:2) = 0;
    side = [1 2; 2 3; 3 1];
    side10 = [1 4 5 2; 2 6 7 3;3 8 9 1 ];
    % Loop over triangles.
    for ee = 1:size(tri,1);
        
        vert = pts(tri(ee,:),:);
        % Generate a straight sided control net for the eeth triangle.
        node{ee} = gen_net(vert);
        
        % Check to see if the current triangle has any sides on a global
        % boundary. If it does, go in an replaces the CPs on that side with
        % the CPs generated from extraction.
        
        for ff = 1:length(edge)
            for ss = 1:3
                if tri(ee,side(ss,:)) == edge(ff,:)
                    node{ee}(side10(ss,:),:) = bPts{ff};
                    bflag(ff,1:2) = [ee,ss];
                elseif tri(ee,side(ss,:)) == fliplr(edge(ff,:))
                    node{ee}(side10(ss,:),:) = flipud(bPts{ff});
                    bflag(ff,1:2) = [ee,ss];                    
                end
            end
        end
    end
end


% generate the IEN and NODE arrays
[IEN,NODE] = gen_arrays(node);

return