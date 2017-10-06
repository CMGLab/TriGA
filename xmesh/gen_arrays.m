function [NODE,IEN] = gen_arrays(node)
%---------------------------------gen_arrays-----------------------------------%
% GEN_ARRAYS generates the NODE and element connectivity (IEN) arrays from
% a cell array of local element nodes.
%
% INPUT:
% node: a 1xnel (number_of_elements) cell array containing the local control
%       points for each element in the mesh.
%
% OUTPUT:
% NODE: a number_of_nodesx3 matrix containing the control points of the
%       mesh.
% IEN: a nenxnel (number_of_element_nodes x number_of_elements) matrix
%      containing the element connectivity infomration for the mesh.
%------------------------------------------------------------------------------%

% Initializing variables
nen = size(node{1},1);
nel = numel(node);
ctr = 1;
allNodes = zeros(nen*nel,5);

% Concatenating all of the local nodes into one large array.
for ee = 1:nel
    allNodes(ctr:ctr+nen-1, :) = [node{ee}(:,1:3) , ones(nen,1)*ee, [1:nen]'];
    ctr=ctr+nen;
end

% Number of decimal places to allow.
d = 10;
tol = 10^-d;
allNodes = [round(allNodes(:,1:2)*10^d)/(10^d) allNodes];
allNodes = sortrows(allNodes,[1,2]);
allNodes = allNodes(:,3:end);

% Pre-allocating the NODE and IEN arrays
IEN = zeros(nen,nel);
NODE = zeros(length(allNodes),2);
NODE(1,1:3) = allNodes(1,1:3);
e = allNodes(1,4);
% Local node number
n = allNodes(1,5);
IEN(n,e) = 1;


ctr = 1;
disp('Generating the connectivity and global node matrices');

for i  = 2:length(allNodes)
    % Element number
    e = allNodes(i,4);
    
    % Local node number
    n = allNodes(i,5);
    
    base_i = i;
    prev_i = i - 1;
        while(  prev_i > 0 && abs( allNodes( i,1) - allNodes( prev_i, 1) ) < tol )
            if( abs( allNodes( i,2) - allNodes( prev_i, 2) ) < tol)
                base_i = prev_i;
            end
            prev_i = prev_i - 1;
        end
    
    if base_i == i
        ctr = ctr+1;
        NODE(ctr,:) = allNodes(i,1:3);
        IEN(n,e) = ctr;
    else
        e0 = allNodes( base_i, 4 );
        n0 = allNodes( base_i, 5 );
        IEN(n,e) = IEN(n0,e0);
    end
end

NODE = NODE(1:ctr,:);

return