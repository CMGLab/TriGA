function [EDGE] = IEN2EDGE(IEN)

nel = size(IEN,2);



allside = zeros(3*nel,5);
side = [1 4 5 2; 2 6 7 3; 3 8 9 1];

ctr = 0;
for ee = 1:nel
    for ss = 1:3
      
        ctr = ctr+1;
        allside(ctr,:) = [IEN(side(ss,:),ee)' ee];
        
    end
end

allside = allside(1:ctr,:);
allside = unique(allside,'rows');
allside = [sort(allside(:,1:4),2) allside];
allside = sortrows(allside,1:4);



EDGE(1,:) = [allside(1,5:9), 0];

ctr = 1;

for ii = 2:nel*3

    if all(allside(ii,1:4) == allside(ii-1,1:4))
        EDGE(ctr,6) = allside(ii,9);
        
    else
        ctr = ctr + 1;
        EDGE(ctr,:) = [allside(ii,5:9), 0];
    end
end

return
        
        