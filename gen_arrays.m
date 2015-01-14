function [IEN, NODE] = gen_arrays(node)
%---------------------------------gen_arrays-----------------------------------%
% GEN_ARRAYS generates the element connectivity (IEN) array from a cell array of
% element nodes. It also generates a list of global nodes, NODE.
%------------------------------------------------------------------------------%

nen = size(node{1},1);
nel = numel(node);
ctr = 1;
allNodes = zeros(nen*nel,4);

for i = 1:nel  
    allNodes(ctr:ctr+nen-1, :) = [node{i}(:,1:2) , ones(nen,1)*i, [1:nen]'];
    ctr=ctr+nen;
end



% Number of decimal places to allow.
d = 10;
allNodes = round(allNodes*10^d)/(10^d);


allNodes = sortrows(allNodes,[1,2]);
% Pre-allocating the IEN array

IEN = zeros(nen,nel);
NODE = zeros(length(allNodes),2);
NODE(1,1:2) = allNodes(1,1:2);
e = allNodes(1,3);
% Local node number
n = allNodes(1,4);
IEN(n,e) = 1;


ctr = 1;
disp('Generating the connectivity and global node matrices');

for i  = 2:length(allNodes)
    % Element number
    e = allNodes(i,3);
    
    % Local node number
    n = allNodes(i,4);
    
 if allNodes(i,1)==allNodes(i-1,1) && allNodes(i,2)==allNodes(i-1,2)
    e0 = allNodes(i-1,3);
    n0 = allNodes(i-1,4);
     IEN(n,e) = IEN(n0,e0);
     
 else
     ctr = ctr+1;
     NODE(ctr,:) = allNodes(i,1:2);
     IEN(n,e) = ctr;
 end

end

NODE = NODE(1:ctr,:);
return