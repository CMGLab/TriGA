function [x] = evalNURBS(P,KV,xi)
%---------------------------------evalNURBS-------------------------------------%
% EVALNURBS: This function evaluates a NURBS curve at a given knot
% location(s) and returns the corresponding coordinates in 2_D physical space.

% INPUT:
% P: nx3 matrix representing the control points of the NURBS Curve in 2-d
% KV: 1xm vector representing the knot vector of the NURBS curve
% xi: Nx1 vector of locations at which to evaluate the NURBS curve.

% OUTPUT:
% x: Nx2 matrix representing the physuical coordinates corresponding to
% each xi. 
%------------------------------------------------------------------------------%

% Getting information about the number of control points, polynomial degree
% and knot vector. 
n = size(P,1);
m = length(KV);
p = m-n-1;


k  = size(P,1);
x = zeros(length(xi),2);

% Looping through and evaluating the curve at each xi using the Cox-deBoor
% recursion relation. 
for j = 1:length(xi)
    num = zeros(1,2);
    den = 0;
    for i = 1:k
        num = num+ deBoor(i,p,xi(j),KV)*P(i,3)*P(i,1:2);
        den = den+ deBoor(i,p,xi(j),KV)*P(i,3);
    end
    x(j,:) = num/den;
end

function N = deBoor(i,p,xi,kv)
%-----------------------------------deBoor-------------------------------------%
% This is an implimentation of the Cox-deBoor Recursion formula for
% evaluating the basis functions for generating NURBS curves. 

% Inputs:
% i: The index number of the basis function
% p: The polynomial degree of the basis function.
% xi: The coordinate (in terms of the knot vector) at which to evaluate the
% function.
% kv: The knot vector for the NURBS curve.

% Output:
% N: The value of the function evaluated at xi.


%------------------------------------------------------------------------------%

if p==0
    
    if xi==kv(i) && xi==kv(i+1)
        N = 0;
    elseif xi>=kv(i) && xi<=kv(i+1)
        N = 1;
    else
        N = 0;
    end

else
    if kv(i) ==kv(i+p)
        f = 0;
    else
        f = (xi - kv(i))/(kv(i+p)-kv(i));
    end
    
    if kv(i+1) == kv(i+p+1)
        g = 0;
    else
        g = (kv(i+p+1)-xi)/(kv(i+p+1)-kv(i+1));
    end
    
    N = f*deBoor(i,p-1,xi,kv) + g*deBoor(i+1,p-1,xi,kv);
    
end

return