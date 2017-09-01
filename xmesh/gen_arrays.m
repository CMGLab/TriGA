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
%    if ee == 43943
%        keyboard
%    end
    allNodes(ctr:ctr+nen-1, :) = [node{ee}(:,1:3) , ones(nen,1)*ee, [1:nen]'];
    ctr=ctr+nen;
end

% Number of decimal places to allow.
d = 14*3;

%allNodes = [round(allNodes(:,1:2)*10^d) allNodes];
tol = 1000;
xtol = max(abs(allNodes(:,1)))*tol*eps();
ytol = max(abs(allNodes(:,2)))*tol*eps();
%[xy_round] = int64(round(allNodes(:,1:2)*2^d));
%[xy_sorted,permutation] = sortrows([xy_round],[1,2]);
[~, allNodes] = sort_perm_matrix( allNodes, xtol, ytol);
%allNodes = allNodes(permutation,:);

%allNodes = sortrows(allNodes,[1,2]);
%allNodes = allNodes(:,3:end);

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
    
 %if round(allNodes(i,1)*2^d) == round(allNodes(i-1,1)*2^d) && ...
 %           round(allNodes(i,2)*2^d) == round(allNodes(i-1,2)*2^d)
        
 if abs(allNodes(i,1) - allNodes(i-1,1)) < xtol && abs(allNodes(i,2) - allNodes(i-1,2)) < ytol
     
       %keyboard
        
    e0 = allNodes(i-1,4);
    n0 = allNodes(i-1,5);
     IEN(n,e) = IEN(n0,e0);
     if single(allNodes(i,3))~=1.0
         NODE(ctr,3) = allNodes(i,3);
     else
         NODE(ctr,3) = 1.0;
     end
 else
     ctr = ctr+1;
     if single(allNodes(i,3))~=1.0
         NODE(ctr,3) = allNodes(i,3);
     else
         NODE(ctr,3) = 1.0;
     end
     NODE(ctr,1:2) = allNodes(i,1:2);
     IEN(n,e) = ctr;
 end

end

NODE = NODE(1:ctr,:);
% 
% for i  = 2:ctr
%     for j = 1:i-1
%         tt = true;
%         for k = 1:2
%             if abs(NODE(i,k) - NODE(j,k)) > max( abs(NODE(i,k)), abs(NODE(j,k)))*10*eps()
%                 tt=false;     
%             end
%         end
%         if tt==true
%             keyboard
%         end
%     end
% end

return