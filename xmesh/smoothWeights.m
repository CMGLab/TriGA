function [NODE] = smoothWeights(NODE, IEN, BFLAG)
% -----------------------------------------------------------------------------%
% SMOOTHWEIGHTS solves a Laplace problem over the mesh to obtain a smoother
% distribution of weights.

% Input:
% NODE: The global node list of the current mesh. 
% IEN: The element connectivity of the current mesh. 
% BFLAG: A list of global node numbers that are on the meshi boundary. 


% Output: 
% NODE: The updates global node list with smoothed weights. 

% -----------------------------------------------------------------------------%

% Rational or Bernstein basis functions
rational = false;

% Set  unit conductivity:
D = eye(2);

% Generating a list of nodes that fall on diriclet boundary conditions.
BFLAG = BFLAG(BFLAG(:,1)~=0,:);
delem = BFLAG(BFLAG(:,4)==6 | BFLAG(:,4)==7,:);
tempBounds = [];
side10  = [1 4 5 2; 2 6 7 3; 3 8 9 1];
for dd = 1:length(delem)
    ee = delem(dd,1);
    ss = delem(dd,2);
    addtempBound = IEN(side10(ss,:),ee);
    tempBounds =[tempBounds;addtempBound]; %#ok<AGROW>
end
    tempBounds = unique(tempBounds);
    boundTemp = NODE(tempBounds,3);


% Setting Glabal Mesh Information.
nel = size(IEN,2);
nen = 10;
nNodes = length(NODE);


% Load the data for the quadrature rule. 
[qPts, qPtsb, W, Wb] = quadData(28);
nQuad   = length(W);
nbQuad  = length(Wb);

% 2D FEA Sovler
% Data initalization. Set global stiffness and forcing matrices to 0.
K(nNodes,nNodes) = sparse(0);
F(nNodes,1) = sparse(0);


% Loop through the elements
for ee =1:nel
    % Display progress in 1% completion increments
    dispFreq = round(nel/100);
    if mod(ee,dispFreq) == 0
        clc
        fprintf('smoothWeights is %3.0f percent complete with the smoothing process\n',ee/nel*100)
    end
    
    % initialize the local stiffness  matrix
    k = zeros(nen);
    
    % Define the local node.
    node = NODE(IEN(:,ee),:);
    
    % Conductivity and heat generation contributions to K and F, respectivly
    % Loop though the quadtrature points.
    for nq = 1:nQuad
        % Call the finite element subroutine to evaluate the basis functions
        % and their derivatives at the current quad point.
        [R, dR_dx,x, J_det]  = tri10(qPts(nq,1), qPts(nq,2),node,rational);
        J_det = abs(J_det);
        
        
        % Add the contribution of the current quad point to the local element
        % stiffness and forcing matrices.
        k = k + W(nq)*dR_dx*D*dR_dx'*J_det;
    end
    
    % Assemble k and f to the glabal matrices, K and F.
    for b = 1:nen
        for a = 1:b
            K(IEN(a,ee),IEN(b,ee)) = K(IEN(a,ee),IEN(b,ee)) + k(a,b);
        end
    end
end


% Fill in the bottom part of the the K matrix .
K = K + K' - K.*speye(size(K));

% ----------Go Back and take care of dirichlet boundary conditions-------------%
disp(' Taking care of diriclet BCs...')

[Ki,Kj,~] = find(K);
for ii = 1:length(tempBounds)
    col = tempBounds(ii);
    row = Ki(Kj == col);
    for j = 1:length(row)
        if row(j)~=col
            F(row(j)) = F(row(j))-K(row(j),col)*boundTemp(ii);
        end
    end
    
    K(col,:) = zeros(1,length(K));
    K(:,col) = zeros(length(K),1);
    
    K(tempBounds(ii),tempBounds(ii)) = 1;
    F(tempBounds(ii)) = boundTemp(ii);
    
end

disp('Solving the System...')
% Solve the system.
temp = K\F;
NODE(:,3) = full(temp);
return
