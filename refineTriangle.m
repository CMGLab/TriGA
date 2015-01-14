function node4 = refineTriangle(node) 

% Generate the locations at which to evaluate tri10 in parameric space.

node3D(:,1) = node(:,1).*node(:,3);
node3D(:,2) = node(:,2).*node(:,3);
node3D(:,3) = node(:,3);

Xi = [0 0;...
    1 0;...
    0   1;...
    1/3 0;...
    2/3 0;...
    2/3 1/3;...
    1/3 2/3;...
    0   2/3;...
    0   1/3;...
    1/3 1/3];

Xi1 = Xi*2;
Xi2 = Xi*2 - [ones(10,1),zeros(10,1)];
Xi3 = Xi*2 - [zeros(10,1),ones(10,1)];
Xi4 = 2*[1 1;...
        0 1;...
        1 0;...
        2/3 1;...
        1/3 1;...
        1/3 2/3;...
        2/3 1/3;...
        1 1/3;...
        1 2/3;...
        2/3 2/3]- [ones(10,1),ones(10,1)];
    
    
for i = 1:10
    Rhat1(i,:) = [tri10(Xi1(i,1),Xi1(i,2),node)]';
    Rhat2(i,:) = [tri10(Xi2(i,1),Xi2(i,2),node)]';
    Rhat3(i,:) = [tri10(Xi3(i,1),Xi3(i,2),node)]';
    Rhat4(i,:) = [tri10(Xi4(i,1),Xi4(i,2),node)]';

end


for i = 1:10
    R(i,:) = [tri10(Xi(i,1),Xi(i,2),node3D)]';
end



x1 = inv(Rhat1)*R*node3D;
x2 = inv(Rhat2)*R*node3D;
x3 = inv(Rhat3)*R*node3D;
x4 = inv(Rhat4)*R*node3D;



x1(:,1) = x1(:,1)./x1(:,3);
x1(:,2) = x1(:,2)./x1(:,3);
x2(:,1) = x2(:,1)./x2(:,3);
x2(:,2) = x2(:,2)./x2(:,3);
x3(:,1) = x3(:,1)./x3(:,3);
x3(:,2) = x3(:,2)./x3(:,3);
x4(:,1) = x4(:,1)./x4(:,3);
x4(:,2) = x4(:,2)./x4(:,3);


node4{1} = x1;
node4{2} = x2;
node4{3} = x3;
node4{4} = x4;

return


