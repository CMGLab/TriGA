function [L2, H1, elemL2,elemL2Rel] =calcNorm(filename,func,gradFunc)
% ---------------------------------------------------------------------------- %
% CALCNORM is a function for calculating the L2 and H1 norms of a
% finite element solution that has a known analytical solution. It is
% mainly used for debugging, testing and evaluating the convergence of
% heat2d. 
%
% INPUT:
% filename: The filename of the gambit  neutral file containing the mesh
% information and the temperature solution field. 
%
% func: a MATLAB function handle contianing the analytical solution as a
% function of cartesian coordinates, (x,y). 
%
% gradFunc: a MATLAB function handle contianing the gradient of the
% analytical solution, again as a function of (x,y). 
%
% OUTPUT: 
%
% L2: The L2 norm of the approximate solution. 
%
% H1: THe H1 norm of the approximate solution. 
%
% elemL2: a nelx1 array containing the L2 norm over each element. This is
% mainly used for debugging purposes. 
%
% elemL2Rel:  a nelx1 array containing the relative L2 norm over each 
% element. The relative error (as opposed to absolute) is calculated and it
% is normalized by element size. This is mainly used for debugging purposes. 
% ---------------------------------------------------------------------------- %

% Reading in the problem data. 
[NODE,IEN,~,temp] = gambitFileIn(filename);

% Initializing variables. 
nel = size(IEN,2); 
nen = size(IEN,1);
elemL2 = zeros(nel,1);
elemL2Rel = zeros(nel,1);
L2 = 0;
H1 = 0;
[qPts, ~, W, ~]  = quadData(28);
nQuad = length(W);
% Generating Lookup Tables for the basis.

[N,dN_du] = evaluateBasis(qPts);
for ee = 1:nel
     % Generating the local node array. 
    node = NODE(IEN(:,ee),:);
    
    % Outputting solver progress. 
    if mod(ee,10) == 0
    clc
        fprintf('calcNorm is %3.0f percent complete\r',ee/nel*100)
    end
    
    L2square = 0;
    H1square = 0;
    sumU = 0;
    for qq = 1:nQuad
        % Find global x location of current quad point
        [R,dR_dx,x,J_det] = tri10fast(node,N(:,qq),dN_du(:,:,qq));
        
        % Evaluate the explicit function at the current quad point.
        u = func(x(1),x(2));
        
        % Evaluate the explicit gradient at the current quad point.
        gradU = gradFunc(x(1),x(2));
        
        % Calculate uh at the current quadrature point.
        uh = 0;
        for i = 1:nen
            uh = uh + R(i)*temp(IEN(i,ee));
        end
        
        % Calculate gradUh at the current quadrature point.
        gradUh = 0;
        for  i  = 1:nen
            gradUh = gradUh + dR_dx(i,:)*temp(IEN(i,ee));
        end

        L2square = L2square + (u-uh)^2*W(qq)/2*J_det;
        H1square = H1square + sum((gradU-gradUh).^2)*W(qq)/2*J_det;
        
        sumU = sumU + u;

    end
    elemL2(ee) = L2square;
    elemL2Rel(ee) = L2square/J_det/(sumU/nQuad);


    L2 = L2 + L2square;
    H1 = H1 + H1square;
end

L2 = sqrt(L2);
H1 = sqrt(H1);
return