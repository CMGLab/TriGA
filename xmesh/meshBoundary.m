function [bNode, BFLAG] = meshBoundary(P, KV,bflag,kvloc)
%---------------------------------meshBoundary---------------------------------%
% This funciton meshes the boundaries of the domain. It takes as an input a
% NURBS curve ( in the form of its control points and knot vector) and a
% user defined number of elements to mesh the edge with. It outputs a series
% of control nets and nodal points corresponding to one "layer" of
% triangular elements along the boundary. It also outputs a list of nodes
% that fall on the outer boundary jsut created (bNode) and a connectivity
% array of these nodes(bEdge). These are for passing to a third party mesh
% generator to mesh the interior volume that's away from the edges.

% This function is effectively Bezier extraction


% Inputs:
% P = control points of the nurbs curve
% KV = knot vector of the nurbs curve.

% Output
% B: a 1 x 2*nel cell of control nets

% node:

ctr = 1;
nCurves = numel(P);
for cc = 1:nCurves
    kvtmp = kvloc(kvloc>=cc-1 & kvloc<=cc)-(cc-1);
    [bNodeTemp,bflagTemp] = extractCurves(P{cc}, KV{cc},kvtmp);   
    bNode(ctr:ctr+numel(bNodeTemp)-1) = bNodeTemp;

    for ii = 1:length(bflagTemp)
        BFLAG(ctr+ii-1) = bflag{cc}(bflagTemp(ii),2);
    end
    
    ctr = ctr + numel(bNodeTemp);
end


return

function [node, bflag] = extractCurves(P,KV,kvloc)
%------------------------------------------------------------------------------%
% Calculate the degree of the inputed Nurbs curve:
% number of knots:
m = length(KV);
% number of control points:
n = size(P,1);
% From the relation: m = n+p+1: -->
p = m-n-1;

% Normalize the knot vector:
KV = KV/KV(end);

% Save a copy of the original knot vector
KV0 = KV;

% Initialize variables:
n_el = length(unique(kvloc)) - 1;
bEdge = zeros(n_el,2);

a = zeros(length(P),1);
Q = zeros(size(P));


%                      BUILDING TARGET KNOT VECTORS
%-------------------------------------------------------------------------%
% Based on the location that we want knots along the final bounary (kvloc)
% build the final knot vector.
if p ==1
    % Calculate the number of control points along the curve:
    nKV = (length(kvloc)-2) + 4;
    
    % Initialize the desired final knot vector:
    KVF = zeros(1,nKV);
    
    % The knot at 0 will have multiplicity 3, so don't change the first
    % three entries of KVF. The knot at the end will also have multiplicity
    % three, change the last three entries to 1.
    KVF(nKV-p:nKV) = [1 1];
    
    % The middle knots will all have multiplicity of 2, so loop through and
    % assign the rest of the knots
    ctr = 3;
    for kk = 2:length(kvloc)-1
        KVF(ctr)   = kvloc(kk);
        ctr = ctr+1;
    end    
elseif p == 2
    % Calculate the number of control points along the curve:
    nKV = (length(kvloc)-2)*2 + 6;
    
    % Initialize the desired final knot vector:
    KVF = zeros(1,nKV);
    
    % The knot at 0 will have multiplicity 3, so don't change the first
    % three entries of KVF. The knot at the end will also have multiplicity
    % three, change the last three entries to 1.
    KVF(nKV-p:nKV) = [1 1 1];
    
    % The middle knots will all have multiplicity of 2, so loop through and
    % assign the rest of the knots
    ctr = 4;
    for kk = 2:length(kvloc)-1
        KVF(ctr)   = kvloc(kk);
        KVF(ctr+1) = kvloc(kk);
        ctr = ctr+2;
    end
    
elseif p == 3
    % Calculate the number of control points along the curve:
    nKV = (length(kvloc)-2)*3 + 8;
    
    % Initialize the desired final knot vector:
    KVF = zeros(1,nKV);
    
    % The first knot will have multiplicity of 4, so don't change the first
    % four entries. The last knot will also have multiplicity of 4, so
    % assign the last four entries to 1
    KVF(nKV-p:nKV) = [1 1 1 1];
    
    % The middle knots will all have multiplicity of 3, so loop through and
    % assign the rest of the knots
    ctr = 5;
    for kk = 2:length(kvloc)-1
        KVF(ctr)   = kvloc(kk);
        KVF(ctr+1) = kvloc(kk);
        KVF(ctr+2) = kvloc(kk);
        ctr = ctr+3;
    end
    
else
    display('This function only supports NURBS curves of degree 2 or 3')
    return
end

%                            KNOT INSTERTION
%-------------------------------------------------------------------------%
% Convert the NURBS curve to a 4D b-spline by multipling all the control
% points by their respective weights.
P(:,1) = P(:,1).*P(:,3);
P(:,2) = P(:,2).*P(:,3);

% Knot Insertion:
% We know what we want our final knot vector to look like, so we'll do knot
% insertion until our KV matches our KVF.
for i = 1:length(KVF)
    % Check to see if the current knot in KV matches the current knot in
    % KVF. If it does, great, move on to the next knot, if not, calculate
    % new control points and insert KVF(i) into KV(i).
    if KVF(i) == KV(i)
        continue
    else
        
        % Correct for the algorithm indexing by zero
        ki = i-1;
        
        % Calculate the new control points
        % Loop through indexes from k-p+1 to k and calculate the new control
        % points
        for j = ki-p+1:ki
            a(j) = (KVF(i)-KV(j))/(KV(j+p)-KV(j));
            Q(j,:) = (1-a(j))*P(j-1,:) + a(j)*P(j,:);
        end
        
        % Move control points below ki-1 down by 1 space to
        % make room for the p new points;
        P(ki+1:length(P)+1,:) = P(ki:length(P),:);
        % Insert the new control points
        P(ki-p+1:ki,:) = Q(ki-p+1:ki,:);
        
        % move the knots to the right of the knot to be inserted over one
        % index
        KV(i+1:end+1) = KV(i:end);
        
        % Insert the knot to the knot vector.
        KV(i) = KVF(i);
    end
end


%                          DEGREE ELEVATION
%-------------------------------------------------------------------------%
% If the inputted curve was quadratic, do degree elevation to get a cubic curve.
if p==1
    for pp = 1:2
    clear Q
    ctr1 = 1;
    ctr2 = 1;
    for e = 1:n_el
        [Q(ctr1:ctr1+p+1,:),pe] = elevateDegree(P(ctr2:ctr2+p,:),p);
        ctr1 = ctr1+p+1;
        ctr2 = ctr2+p;
    end
    
    
    P = Q;
    p  = pe;
    
    kvloc = unique(KV);
    nKV = length(P)+p+1;
    % Go in and change the knot vector.
    % Initialize the desired final knot vector:
    KVe = zeros(1,nKV);
    
    % The first knot will have multiplicity of 4, so don't change the first
    % four entries. The last knot will also have multiplicity of 4, so
    % assign the last four entries to 1
    KVe(nKV-p:nKV) = ones(1,p+1);
    
    % The middle knots will all have multiplicity of 3, so loop through and
    % assign the rest of the knots
    ctr = p+1;
    for kk = 2:length(kvloc)-1
        KVe(ctr:ctr+p-1)   = kvloc(kk);
        ctr = ctr+p;
    end
    KV = KVe;
    end
elseif p==2
    clear Q
    ctr1 = 1;
    ctr2 = 1;
    for e = 1:n_el
        [Q(ctr1:ctr1+3,:),pe] = elevateDegree(P(ctr2:ctr2+2,:),p);
        ctr1 = ctr1+3;
        ctr2 = ctr2+2;
    end
    
    
    P = Q;
    p  = pe;
    
    kvloc = unique(KV);
    nKV = length(P)+p+1;
    % Go in and change the knot vector.
    % Initialize the desired final knot vector:
    KVe = zeros(1,nKV);
    
    % The first knot will have multiplicity of 4, so don't change the first
    % four entries. The last knot will also have multiplicity of 4, so
    % assign the last four entries to 1
    KVe(nKV-p:nKV) = [1 1 1 1];
    
    % The middle knots will all have multiplicity of 3, so loop through and
    % assign the rest of the knots
    ctr = 5;
    for kk = 2:length(kvloc)-1
        KVe(ctr)   = kvloc(kk);
        KVe(ctr+1) = kvloc(kk);
        KVe(ctr+2) = kvloc(kk);
        ctr = ctr+3;
    end
    KV = KVe;
end
%-------------------------------------------------------------------------%
%                         END DEGREE ELEVATION

% Now that the knot insertion is done, renormalize all the control points
% by their respective weights to get a NURBS description again.
P(:,1) = P(:,1)./P(:,3);
P(:,2) = P(:,2)./P(:,3);
%-------------------------------------------------------------------------%
%                         END KNOT INSTERTION


d = 14;

bb = 1:3:length(P);
bNode = P(bb(1:end-1),1:2);
for ee = 1:n_el
    bEdge(ee,:) = [ee,ee+1];
    node{ee} = round(P(bb(ee):bb(ee+1),:)*10^d)/10^d;
end

bEdge(end,2) = 1;

% Find out which knot span of the original knot vector each element lies
% on
KV0 = unique(KV0);
KVF = unique(KV);
bb = zeros(n_el,1);
for kk = 2:length(KV0)
    bb(kk) = sum(KVF>KV0(kk-1) & KVF<KV0(kk))+1;
    bsum = cumsum(bb);
    bflag(bsum(kk-1)+1:bsum(kk),1) = ones(bb(kk),1)*(kk-1);
end


return

function [Pe,pe] = elevateDegree(P,p)

% This function performs degree elevation on a NURBS curve of degree p. 
% It takes as input the control points of the existing nurbs curve, P and its
% knot vector and outputs the control points and knot vector if a degree p+1
% NURBS curve. 


n = size(P,1);

ne = n+1;
pe = p+1;

Pe(1,:) = P(1,:);
Pe(ne,:) = P(n,:);

for i = 2:n
    Pe(i,:) = (i-1)/(p+1)*P(i-1,:) + (1 - (i-1)/(p+1))*P(i,:);
end

return