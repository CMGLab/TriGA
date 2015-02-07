function [node] = gen_net(vert)
% ----------------------------------gen_net------------------------------------%
% Generates the control net for the current triangular Bezier patch. Reads in
% the verticies of a triangle, and linearly interpolates the triangle to make a
% 10 control point bezier element. 

% Input:
% vert: a matrix in the form [x1 y1; x2 y2; x3 y3] containing the verticies of a
%  linear triangle, ordered counter clockwise from the bottom left vertex.

%     v3
%     | \
%     |  \
%     |   \
%     |    \
%     v1----v2

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
%     1--4--5--2
%
%------------------------------------------------------------------------------%

nen = 10;
idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];

node(:,1:2) = [ vert(1,1:2);...
         vert(2,1:2);...
         vert(3,1:2);...
         vert(1,1:2) + (vert(2,1:2)-vert(1,1:2))/3;...
         vert(1,1:2) + (vert(2,1:2)-vert(1,1:2))*2/3;...
         vert(2,1:2) + (vert(3,1:2)-vert(2,1:2))/3;...
         vert(2,1:2) + (vert(3,1:2)-vert(2,1:2))*2/3;...
         vert(3,1:2) + (vert(1,1:2)-vert(3,1:2))/3;...
         vert(3,1:2) + (vert(1,1:2)-vert(3,1:2))*2/3];
            
node(10,1:2) = node(9,1:2) + (node(6,1:2)-node(9,1:2))/2;
if size(vert,2) ==3
node(4:10,3) = ones(7,1);
node(1:3,3) = vert(1:3,3);
else
    node(:,3) = ones(10,1);
end

return