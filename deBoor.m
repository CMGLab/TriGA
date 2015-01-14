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