function [R, dR_dx,x, J_det] = tri10fast(node,N,dN_du,varargin)

nen = length(N);
if nargin == 3
    rational = true;
elseif nargin == 4
    rational = varargin{1};
end
%-------------------------------------------------------------------------%
% TRI10 is the finite element subroutine for a 10 node triangular element.
% It uses either Rational Beziers or Bernstein Polynomials as basis
% functions over the unit triangle . 
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
%     3
%     |\
%     | \
%     8  7
%     |   \
%     |    \
%     9 10  6
%     |      \
%     |       \
%     1--4--5--2
%
% rational: (optional) If rational == true, rational bezier basis functions
% will be used. If rational == false, Bernstein polynomial will be used. If
% no value is specified, tri10 defaults to Rational Beziers. 
%
% OUTPUT: 
% R: A 10x1 array containing the basis functions evaluated at [xi,eta].
%
% dR_dx: A 10x2 array containing the basis function derivatives 
% [dR_dx, dR_dy] evaluated at [xi,eta]
%
% x: The [x,y] location in physical space corresponding to [xi,eta].
%
% J_det: The Jacobian determinant of the mapping from parametric space to
% physical space.
%------------------------------------------------------------------------------%
R = zeros(nen,1);
dR_du = zeros(nen,3);
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
    dR_du = dN_du;
end


num = 0;
den = 0;
for i = 1:nen
    num = num +N(i)*node(i,1:2)*node(i,3);
    den = den +N(i)*node(i,3);
    
end
x = num/den;

% Chain rule to find the derivative with respect to cartesian isoparametric
% coordinates.
du_dxi = [-1 -1; 1 0;0 1];
dR_dxi = dR_du*du_dxi;
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

dx_dxi = (gp.*h-g.*hp)./(h.^2);

% Calculating the shape function derivatives and the Jacobian determinate.
dR_dx =dR_dxi*inv(dx_dxi); %#ok<MINV>
J_det = det(dx_dxi);
return


