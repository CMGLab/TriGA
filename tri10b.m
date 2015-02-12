function [R, J_detb] = tri10b(xi,eta,node,s,varargin)
%-------------------------------------------------------------------------%
% TRI10B is the finite element subroutine for the edges of a 10 node 
% triangular element. It uses either Rational Beziers or Bernstein 
% Polynomials as basis functions over the unit triangle. 
%
% INPUT:
% xi: The xi location (in parametric space) at which to evaluate the basis
% fucntions.
%
% eta: The eta location (in parametric space) at which to evaluate the basis
% fucntions.
%
% node: A 10x3 array that represents the control net for a 10-node Bezier
% triangle. The first two columns of node give the x and y coordinates of 
% the nodes in physical space, and the last column gives their 
% corresponding weights. Node ordering follows the convention below:
%
%           3
%           |\
%           | \
%           8  7
%           |   \
%   side 3  |    \   side 2
%           9 10  6
%           |      \
%           |       \
%           1--4--5--2
%
%             side 1
%
% s: The side of the triangle on which to evaluate the Basis finctions.
% The side numbering convention is shown above.
%
% rational: (optional) If rational == true, rational bezier basis functions
% will be used. If rational == false, Bernstein polynomial will be used. If
% no value is specified, tri10 defaults to Rational Beziers. 
%
% OUTPUT: 
% R: A 10x1 array containing the basis functions evaluated at [xi,eta].
%
% J_det: The Jacobian determinant of the mapping from parametric space to
% physical space.
%------------------------------------------------------------------------------%

if nargin == 4
    rational = true;
elseif nargin == 5
    rational = varargin{1};
end

% Element parameters.
n = 3;
nen = 10;
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

% Initializing variables
N = zeros(nen,1);
R = zeros(nen,1);
dN_du = zeros(nen,3);
dR_du = zeros(nen,3);

du_dxi = [-1 -1; 1 0;0 1];

% Index and tuples. The index is the location in the control net (row,col) of
% the ith control point. Tuples is the index in barycentric coordinates of the
% ith control point.
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
    N(nn) = factorial(n)/...
        (factorial(i)*factorial(j)*factorial(k))*u^i*v^j*w^k;
    
    if i-1 < 0
        dN_du(nn,1) =0;
    else
        dN_du(nn,1) = n * factorial(n-1)/...
            (factorial(i-1)*factorial(j)*factorial(k))*u^(i-1)*v^j*w^k;
    end
    
    if j-1 < 0
        dN_du(nn,2) = 0;
    else
        dN_du(nn,2) = n * factorial(n-1)/...
            (factorial(i)*factorial(j-1)*factorial(k))*u^(i)*v^(j-1)*w^k;
    end
    
    if k-1 < 0;
        dN_du(nn,3) = 0;
    else
        dN_du(nn,3) = n * factorial(n-1)/...
            (factorial(i)*factorial(j)*factorial(k-1))*u^i*v^j*w^(k-1);
    end
    
end

if rational
    den  = N'*node(:,3);
    for nn = 1:nen
        R(nn) = N(nn)*node(nn,3) / den;
        
        dR_du(nn,1) =  (dN_du(nn,1)*node(nn,3)*den - ...
            dN_du(:,1)'*node(:,3)*N(nn)*node(nn,3))/den^2;
        dR_du(nn,2) =  (dN_du(nn,2)*node(nn,3)*den - ...
            dN_du(:,2)'*node(:,3)*N(nn)*node(nn,3))/den^2;
        dR_du(nn,3) =  (dN_du(nn,3)*node(nn,3)*den - ...
            dN_du(:,3)'*node(:,3)*N(nn)*node(nn,3))/den^2;
    end
else
    R = N;
end

% Chain rule to find the derivative with respect to cartesian isoparametric
% coordinates.
dN_dxi = dN_du*du_dxi;

% Calculating the mapping from isoparametric space to physical space.
g = [0 0;0 0];
gp = [0 0; 0 0 ];
h = [0 0  ;0 0];
hp = [0 0; 0 0 ];

for row  = 1:2
    for col = 1:2
        for nn = 1:nen
            g(row,col) = g(row,col) + N(nn)*node(nn,row)*node(nn,3);
            h(row,col) = h(row,col) + N(nn)*node(nn,3);
            gp(row,col) = gp(row,col) + dN_dxi(nn,col)*node(nn,row)*node(nn,3);
            hp(row,col) = hp(row,col) + dN_dxi(nn,col)*node(nn,3);
        end
    end
end
dx_dxi = (gp.*h-g.*hp)./h.^2;

if s == 1
    J_detb = sqrt(sum(dx_dxi(:,1).^2));
elseif s == 2
    J_detb = sqrt((dx_dxi(1,1)-dx_dxi(1,2))^2 + (dx_dxi(2,1)-dx_dxi(2,2))^2);    
elseif s == 3
    J_detb = sqrt(sum(dx_dxi(:,2).^2));
end
return


