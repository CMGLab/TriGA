function [L2, H1, elemL2,elemL2Rel] =calcNorm(inFile,func,gradFunc)
% ---------------------------------------------------------------------------- %


% ---------------------------------------------------------------------------- %

load(inFile)
load quadData28
nel = numel(node); 
nen = size(node{1},1);
elemL2 = zeros(nel,1);
elemL2Rel = zeros(nel,1);
L2 = 0;
H1 = 0;
nQuad = length(W);
for e = 1:nel
            
    if mod(e,10) == 0
    clc
        fprintf('calcNorm is %3.0f percent complete\r',e/nel*100)
    end
    
    L2square = 0;
    H1square = 0;
    sumU = 0;
    for q = 1:nQuad
        % Find global x location of current quad point
        [R,dR_dx,detJ] = tri10(qPts(q,1),qPts(q,2),node{e});
        
        num = 0;
        den = 0;        
        for i = 1:nen
            num = num +R(i)*node{e}(i,1:2)*node{e}(i,3);
            den = den +R(i)*node{e}(i,3);

        end
        x = num/den;
        
        % Evaluate the explicit function at the current quad point.
        u = func(x(1),x(2));
        
        % Evaluate the explicit gradient at the current quad point.
        gradU = gradFunc(x(1),x(2));
        
        % Calculate uh at the current quadrature point.
        uh = 0;
        for i = 1:nen
            uh = uh + R(i)*temp(IEN(i,e));
        end
        
        % Calculate gradUh at the current quadrature point.
        gradUh = 0;
        for  i  = 1:nen
            gradUh = gradUh + dR_dx(i,:)*temp(IEN(i,e));
        end

        L2square = L2square + (u-uh)^2*W(q)/2*detJ;
        H1square = H1square + sum((gradU-gradUh).^2)*W(q)/2*detJ;
        
        sumU = sumU + u;

    end
    elemL2(e) = L2square;
    elemL2Rel(e) = L2square*detJ/(sumU/nQuad);


    L2 = L2 + L2square;
    H1 = H1 + H1square;
end

L2 = sqrt(L2);
H1 = sqrt(H1);
save(inFile,'L2','H1','elemL2','elemL2Rel','-append')
return