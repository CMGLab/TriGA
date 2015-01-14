function [] = showResults(fileName)
% -----------------------------------------------------------------------------%
% SHOWRESULTS reads in the data from fea2d stored in 'fileName.mat' and prints a
% heatmap of the results.
% -----------------------------------------------------------------------------%

load(fileName)

idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];
for ee  = 1:numel(node)
    T1 = node{ee};
    for i = 1:10
        B{ee}(idx(i,1),idx(i,2),1) = T1(i,1);
        B{ee}(idx(i,1),idx(i,2),2) = T1(i,2);
    end
    
end


figure
hold on
nen  = size(node{1},1);
if nen ==10
    for i = 1:numel(node)
        clc
        fprintf('Plotting element %4.0f of %4.0f \n',[i,numel(node)])
        gen_patch(B{i},temp(IEN(:,i)),1/30)
    end
    
elseif nen == 3
    for i = 1:numel(node)
        clc
        fprintf('Plotting element %4.0f of %4.0f \n',[i,numel(node)])
        gen_patchLin(node{i},temp(IEN(:,i)))
    end
    
end
    axis([min(NODE(:,1)) max(NODE(:,1)) min(NODE(:,2)) max(NODE(:,2))])
    colorbar
    
    return
    
