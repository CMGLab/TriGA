function [NODE,EDGE,KVLOC,BFLAG] = NURBS2poly(P,KV,thresh)

%-------------------------------------------------------------------------%
%NURBS2POLY This function approximates any number of closed NURBS curves by
% an equal number of polygons.

% INPUTS:
% P: 1XnCurves cell array containing the control points for each curve
% KV: 1xnCurves cell array contianing the knot vectors for each curve.
% thresh: A threshold for how closely the polygon must match the inputed
% spline geometry. This is set at a default of 1.01 (the lengths of the
% polygon and spline must be within 1% of each other).

% OUTPUTS:
% node: nPtsx2 array containing the coordinates of the polygons.
% edge: nEdgesx2 array containing the connectivity of the polygons.
%-------------------------------------------------------------------------%


% Start by looping through the NURBS curves.
nCurves = numel(P);

CTR = 1;
KVLOC = [];
for cc = 1:nCurves
    clear node kvloc addkv addnode edge
    % Start by forming a polygon by linearly interpolating all of the knot
    % spans and their midpoints.
    
    
    % Building the location knot vector. This array just contains the
    % locations of the knots, but no multiplicity information. NURBS2poly
    % adds a knot to every knot span by default. This has been shown to
    % produce better meshes in general.
    kvloc = unique(KV{cc});
    addkv = zeros(1,length(kvloc)-1);
    for kk = 1:length(kvloc)-1;
        addkv(kk) = (kvloc(kk) + kvloc(kk+1))/2;
    end
    kvloc = sort([kvloc,addkv]);
    
    
    % Loop through the knot spans and evaluate the NURBS curve at the given
    % knot locations. These points become the vertices of the polygon.
    node = zeros(length(kvloc)-1,2);
    edge = zeros(length(kvloc)-1,2);
    curveLength = zeros(length(kvloc)-1,1);
    edgeLength  = zeros(length(kvloc)-1,1);
    for kk = 1:length(kvloc)-1
        node(kk,:)  = evalNURBS(P{cc},KV{cc},kvloc(kk));
        edge(kk,:) = [kk kk+1];
        curveLength(kk) = calcArc(P{cc},KV{cc},[kvloc(kk), kvloc(kk+1)]);
        
        p1 = node(kk,:);
        p2 = evalNURBS(P{cc},KV{cc},kvloc(kk+1));
        edgeLength(kk) = sqrt(sum((p2-p1).^2));
    end
    node(end+1,:) = node(1,:); %#ok<AGROW>
    edge(end)  =1;
    
    % Now check to see which edges we need to split.
    split = curveLength./edgeLength > thresh;
    
    maxit = 10; % Maximum number of iterations to allow. 
    it = 1;
    while true && it<=maxit
        if any(split)
            idx = find(split);
            addkv = zeros(1,length(idx));
            addnode = zeros(length(idx),2);
            % Looping through sides to split.
            
            for ii = 1:length(idx)
                p1 = edge(idx(ii),1);
                p2 = edge(idx(ii),2);
                
                %Inserting new nodes
                addkv(ii) = (kvloc(p1)+kvloc(p2))/2;
                addnode(ii,:) = evalNURBS(P{cc},KV{cc},addkv(ii));
                
                edge(p1+1:end+1,:) = edge(p1:end,:)+1;
                edge(p1,:) = [p1 p1+1];
                
            end
            augM = [[node;addnode], [kvloc(:);addkv(:)]];
            augM = sortrows(augM,3);
            
            node = augM(:,1:2);
            kvloc = augM(:,3);
            edge(end) = 1;
            
        else
            break
        end
        
        % Calculate the new edge and curve lengths and check again to see
        % if any edges need to be split. 
        curveLength = zeros(length(kvloc)-1,1);
        edgeLength  = zeros(length(kvloc)-1,1);
        for kk = 1:length(kvloc)-1
            curveLength(kk) = calcArc(P{cc},KV{cc},[kvloc(kk), kvloc(kk+1)]);
            p1 = node(kk,:);
            p2 = node(kk+1,:);
            edgeLength(kk) = sqrt(sum((p2-p1).^2));
        end
        
        split = curveLength./edgeLength > thresh;
        it = it+1;
    end
    
    if it == 11
        disp('Warning: The maximum number of iterations for NURBS2Poly was exceeded.')
        disp('Please email unmeshable geometries to engvall@colorado.edu')
    end
    
    % Save out the local curve information to the global matrices
    n_el = length(node)-1;
    NODE(CTR:CTR+n_el-1,:) = node(1:end-1,:);
    EDGE(CTR:CTR+n_el-1,:) = edge+CTR-1;
    CTR  = CTR +n_el;
    
    KVLOC= unique([KVLOC kvloc(:)'/kvloc(end)+cc-1]);
    BFLAG = ones(length(kvloc)-1,1);
end


return


function [s] = calcArc(P,KV,kvspan)
% Approximates the arclength of the curve by subdividing into 10 increments and
% suming the linears.
xi  = linspace(kvspan(1),kvspan(2),11);
C = evalNURBS(P,KV,xi);
s = 0;
for i  = 2:length(C);
    dx = C(i,1)-C(i-1,1);
    dy = C(i,2)-C(i-1,2);
    s = s + sqrt(dx^2 + dy^2);
    
end
return