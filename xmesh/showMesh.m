function [] = showMesh(filename)
%-------------------------------------------------------------------------%
% SHOWMESH This function displays the mesh described by the data in
% <filename>.neu. 
%
% INPUT: 
% filename: The name of the gambt neutral file containing the mesh
% information. 
%-------------------------------------------------------------------------%

[NODE,IEN] = gambitFileIn(filename);

figure('Name','Mesh')
hold on
axis equal off
n_el = size(IEN,2);

xi = 0:.1:1;
side10 = [1 4 5 2; 2 6 7 3;3 8 9 1 ];
for ee =  1:n_el
    node = NODE(IEN(:,ee),:);
    for ss = 1:3;
        C = evalNURBS(node(side10(ss,:),:),[0 0 0 0 1 1 1 1],xi);
        plot(C(:,1),C(:,2),'b')
    end
end

return