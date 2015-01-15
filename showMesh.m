function [] = showMesh(NODE,IEN)

% THis function takes as input the cell array of local control points and
% outputs a plot of the mesh. 

figure
hold on
axis equal

n_el = size(IEN,2);
side10 = [1 4 5 2; 2 6 7 3;3 8 9 1 ];
for ee =  1:n_el
    node = NODE(IEN(:,ee),:);
    for ss = 1:3;
        C = gen_curve(node(side10(ss,:),:),[0 0 0 0 1 1 1 1],3);
        plot(C(:,1),C(:,2),'b')
    end
end

return