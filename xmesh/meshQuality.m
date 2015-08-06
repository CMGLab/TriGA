function [J,S] = meshQuality(filename)
%-------------------------------------------------------------------------%
% MESHQUALITY- This function calculates the mesh quality of higher order
% triangular meshes according to two metrics. 

% INPUT:
% filename: filename of the gambit neutral file containing the mesh. 

% OUTPUT: 
% J: A measure of element size variation accross an element.
% S: A measure of element skewness.


%-------------------------------------------------------------------------%

addpath('~/TriGA/');

% Load the mesh information.
[NODE, IEN] = gambitFileIn(filename);
nel = length(IEN);

% Set the quadrature rule.
nq = 28;
[qPts]  = quadData(nq);


% Make lookup tables for the Bernstein basis. 
[N,dN_du] = evaluateBasis(qPts);

% Initialize variables.
dx_dxi = zeros(2,2,nq);
J_det = zeros(nq,1);
J = zeros(nel,1);
S = zeros(nel,1);

% Loop over elements in the mesh. 
for ee = 1:nel 
   node = NODE(IEN(:,ee),:); 
   
   % Loop over quadrature points. 
   
   for qq = 1:nq
        [~, ~,~, J_det(qq),dx_dxi(:,:,qq)] = tri10fast(node,N(:,qq),dN_du(:,:,qq));      
   end
   
   J(ee) = min(abs(J_det))/max(abs(J_det));
   
   S(ee) = max(abs(dx_dxi(:)))/sqrt(max(abs(J_det)));
        
end

return