function [N,dN_du] = evaluateBasis(qPts)

% Element parameters.
n = 3;
nen = 10;
nQPts = size(qPts,1);
N = zeros(nen,nQPts);
dN_du = zeros(nen,3,nQPts);

for qq = 1:nQPts
    xi = qPts(qq,1);
    eta = qPts(qq,2);
    
    % Find the barycentric coordinates of xi and eta
    vert = [0,0;1,0;0,1];
    x1 = vert(1,1);
    x2 = vert(2,1);
    x3 = vert(3,1);
    y1 = vert(1,2);
    y2 = vert(2,2);
    y3 = vert(3,2);
    
    detA =  det([x1,x2,x3; y1,y2,y3;1, 1, 1]);
    detA1 = det([xi,x2,x3;eta,y2,y3;1, 1, 1]);
    detA2 = det([x1,xi,x3;y1,eta,y3;1, 1, 1]);
    detA3 = det([x1,x2,xi;y1,y2,eta;1, 1, 1]);
    
    u = detA1/detA;
    v = detA2/detA;
    w = detA3/detA;
    
    % Tuples is the index in barycentric coordinates of the ith control point.
    tuples  = [ 3 0 0;...
        0 3 0;...
        0 0 3;...
        2 1 0;...
        1 2 0;...
        0 2 1;...
        0 1 2;...
        1 0 2;...
        2 0 1;...
        1 1 1];
    
    % Loop through the control points.
    for nn = 1:nen
        
        i = tuples(nn,1);
        j = tuples(nn,2);
        k = tuples(nn,3);
        
        % From page 141 of Bezier and B-splines. Calculate the ith basis function
        % its derivative with respect to barycentric coordinates.
        N(nn,qq) = factorial(n)/...
            (factorial(i)*factorial(j)*factorial(k))*u^i*v^j*w^k;
        
        if i-1 < 0
            dN_du(nn,1,qq) =0;
        else
            dN_du(nn,1,qq) = n * factorial(n-1)/...
                (factorial(i-1)*factorial(j)*factorial(k))*u^(i-1)*v^j*w^k;
        end
        
        if j-1 < 0
            dN_du(nn,2,qq) = 0;
        else
            dN_du(nn,2,qq) = n * factorial(n-1)/...
                (factorial(i)*factorial(j-1)*factorial(k))*u^(i)*v^(j-1)*w^k;
        end
        
        if k-1 < 0;
            dN_du(nn,3,qq) = 0;
        else
            dN_du(nn,3,qq) = n * factorial(n-1)/...
                (factorial(i)*factorial(j)*factorial(k-1))*u^i*v^j*w^(k-1);
        end
        
    end
end
return
