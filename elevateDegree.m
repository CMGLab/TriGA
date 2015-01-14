function [Pe,pe] = elevateDegree(P,p)

% This function performs degree elevation on a NURBS curve of degree p. 
% It takes as input the control points of the existing nurbs curve, P and its
% knot vector and outputs the control points and knot vector if a degree p+1
% NURBS curve. 

% P(:,1) = P(:,1).*P(:,3);
% P(:,2) = P(:,2).*P(:,3);

n = size(P,1);

ne = n+1;
pe = p+1;

Pe(1,:) = P(1,:);
Pe(ne,:) = P(n,:);

for i = 2:n
    Pe(i,:) = (i-1)/(p+1)*P(i-1,:) + (1 - (i-1)/(p+1))*P(i,:);
end

% Pe(:,1) = Pe(:,1)./Pe(:,3);
% Pe(:,2) = Pe(:,2)./Pe(:,3);


return