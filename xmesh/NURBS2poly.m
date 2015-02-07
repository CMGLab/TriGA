function [NODE,EDGE,KVLOC,BFLAG] = NURBS2poly(P,KV,thresh)

%-------------------------------------------------------------------------%
%NURBS2POLY This function approximates any number of closed NURBS curves by
%an equal amount of polygons.

% INPUTS:
% P: 1XnCurves cell array containing the control points for each curve
% KV: 1xnCurves cell array contianing the knot vectors for each curve.

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
    % Calculating the polynomial degree
    p = length(KV{cc})-length(P{cc}) -1;
    % Start by forming a polygon by linearly interpolating all of the knot
    % spans and their midpoints.
    
    
    % Building the location knot vector,
    kvloc = unique(KV{cc});
    
    for kk = 1:length(kvloc)-1;
        addkv(kk) = (kvloc(kk) + kvloc(kk+1))/2;
    end
    
    kvloc = sort([kvloc,addkv]);
    
    
    ctr = 1;
    
    for kk = 1:length(kvloc)-1
        node(ctr,:)  = evalNURBS(P{cc},KV{cc},kvloc(kk));
        
        edge(ctr,:) = [ctr ctr+1];
        
        cLength(ctr) = calcArc(P{cc},KV{cc},[kvloc(kk), kvloc(kk+1)]);
        
        p1 = node(ctr,:);
        p2 = evalNURBS(P{cc},KV{cc},kvloc(kk+1));
        eLength(ctr) = sqrt(sum((p2-p1).^2));
        ctr = ctr+1;
        
        
    end
    node(end+1,:) = node(1,:);
    edge(end)  =1;
       
    % Now check to see which edges we need to split.
    %     thresh = 1.01;
    
    split = cLength./eLength > thresh;
    clear addkv addnode
    while true
        if any(split)
            idx = find(split);
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
        for kk = 1:length(kvloc)-1
            
            cLength(kk) = calcArc(P{cc},KV{cc},[kvloc(kk), kvloc(kk+1)]);
            
            p1 = node(kk,:);
            p2 = node(kk+1,:);
            eLength(kk) = sqrt(sum((p2-p1).^2));
            
            
            
        end
        
        split = cLength./eLength > thresh;
        
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