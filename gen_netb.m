function [B, node] = gen_netb(bPts,varargin)
% ----------------------------------gen_netb-----------------------------------%
% Generates the control net for the current triangular Bezier patch that is on a
% boundary. Reads in the four control points that describe the curve that lies
% on the boundary. Determines the location of the third vertex fo the triangle
% and generates a net of 10 control points.

% Input:
% vert: a matrix in the form [x1 y1 w1; x2 y2 w2; x3 y3 w3; x4 y4 w4] containing
% the coordinates and weights of the control points describing the curved
% boundary.



% Output: 

%     3
%     |\
%     | \
%     8  7
%     |   \
%     |    \
%     9 10  6
%     |      \
%     |       \
%     1--4--5--2 <------nodes 1,4,5,2 are on the curved boundary.   
%               
%------------------------------------------------------------------------------%
vert(1,:) = bPts(4,1:2);
vert(2,:) = bPts(1,1:2);

if nargin == 1
% Determining the location of the third vertex of the boundary triangle.
    dx = bPts(4,1:2)-bPts(1,1:2);
    mag = sqrt(sum(dx.^2));
    nml = [dx(2), -dx(1)]/mag;
    vert(3,:) = bPts(1,1:2)+dx/2 + nml*mag;

elseif nargin == 2
    
    vert(3,:) = varargin{1}(1:2);
    
else
    disp('Too many input arguments')
    return
end

% Initializing variables.
B = zeros(4,4,3);
B(:,:,3) = ones(4);
nen = 10;
idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];

% Interpolating the sides and generating the control points.
node = [ vert(1,:);...
         vert(2,:);...
         vert(3,:);...
         bPts(3,1:2)
         bPts(2,1:2) 
         vert(2,:) + (vert(3,:)-vert(2,:))/3;...
         vert(2,:) + (vert(3,:)-vert(2,:))*2/3;...
         vert(3,:) + (vert(1,:)-vert(3,:))/3;...
         vert(3,:) + (vert(1,:)-vert(3,:))*2/3];
            
     
% Calculating the location of the center node. Calculated by finding the
% intersection of the line going from node 4 to node 7 with the line going from
% node 5 to node 8.
m1 = (node(7,2)-node(4,2))/(node(7,1)-node(4,1));
m2 = (node(8,2)-node(5,2))/(node(8,1)-node(5,1));
x = (node(5,2)-node(4,2)+m1*node(4,1)-m2*node(5,1))/(m1-m2);
y = m1*(x-node(4,1))+node(4,2);
node(10,:) = [x, y];

% Assigning the weights to the control points. Only 1 4 and 5 2 have weights not
% equal to one:
node(:,3) = ones(10,1);
node([1 4 5 2] ,3) = flipud(bPts(:,3));

for n = 1:nen
    B(idx(n,1),idx(n,2),1:3) = node(n,1:3);
end
return