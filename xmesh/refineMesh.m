function [NODE4,IEN4,BFLAG4] = refineMesh(filename)

% Reading in the data file. 
[NODE,IEN, BFLAG] = gambitFileIn(filename);

% Setting global mesh data. 
nel = size(IEN,2);
ctr = 0;

% Initializing refined mesh variables.
node4 = cell(1,nel*4);
BFLAG4 = zeros(length(BFLAG)*4,4);


% Loop over the elements in the mesh/ 
for ee = 1:nel
    
    % Generate the local node array, and perform uniform refinement using
    % refine Triange.
    node = NODE(IEN(:,ee),:);
    tempNode = refineTriangle(node);
    for nn = 1:4
        node4{ctr+nn} = tempNode{nn};
    end
    
    % Update the boundary condution matrix.
    temp  = BFLAG(BFLAG(:,1)==ee,:);
    for bb = 1:size(temp,1)
    
        if temp(bb,2) == 1
            BFLAG4(ctr+1,:) = [ctr+1 temp(bb,2:4)];
            BFLAG4(ctr+2,:) = [ctr+2 temp(bb,2:4)];
        end
        
        if temp(bb,2) == 2
            BFLAG4(ctr+2,:) = [ctr+2 temp(bb,2:4)];
            BFLAG4(ctr+3,:) = [ctr+3 temp(bb,2:4)];
        end
        
        if temp(bb,2) == 3
            BFLAG4(ctr+1,:) = [ctr+1 temp(bb,2:4)];
            BFLAG4(ctr+3,:) = [ctr+3 temp(bb,2:4)];
        end
    end
    ctr = ctr+4;
end

% Generate the global NODE and IEN arrays from the local node arrays. 
[NODE4,IEN4] = gen_arrays(node4);

% Write out a gambit file
filename = [filename,'ref'];
gambitFileOut(filename,NODE4,IEN4,BFLAG4);

return


function node4 = refineTriangle(node)

nen = size(node,1);

% Transforming the control points to projective space. 
node3D(:,1) = node(:,1).*node(:,3);
node3D(:,2) = node(:,2).*node(:,3);
node3D(:,3) = node(:,3);

% Generate the locations at which to evaluate tri10 in parametric space.
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

Xi1 = Xi*1/2;
Xi2 = Xi*1/2 + 1/2*[ones(10,1),zeros(10,1)];
Xi3 = Xi*1/2 + 1/2*[zeros(10,1),ones(10,1)];
Xi4 = 1/2*[1-Xi(:,1) 1-Xi(:,2)];


Rhat1 =  zeros(nen);
Rhat2 =  zeros(nen);
Rhat3 =  zeros(nen);
Rhat4 =  zeros(nen);
R     =  zeros(nen);
for i = 1:10
    Rhat1(i,:) = tri10(Xi1(i,1),Xi1(i,2))';
    Rhat2(i,:) = tri10(Xi2(i,1),Xi2(i,2))';
    Rhat3(i,:) = tri10(Xi3(i,1),Xi3(i,2))';
    Rhat4(i,:) = tri10(Xi4(i,1),Xi4(i,2))';
    
end


for i = 1:10
    R(i,:) = tri10(Xi(i,1),Xi(i,2))';
end

x1 = R\Rhat1*node3D;
x2 = R\Rhat2*node3D;
x3 = R\Rhat3*node3D;
x4 = R\Rhat4*node3D;


% Converting back from projective space to physical space. 
x1(:,1) = x1(:,1)./x1(:,3);
x1(:,2) = x1(:,2)./x1(:,3);
x2(:,1) = x2(:,1)./x2(:,3);
x2(:,2) = x2(:,2)./x2(:,3);
x3(:,1) = x3(:,1)./x3(:,3);
x3(:,2) = x3(:,2)./x3(:,3);
x4(:,1) = x4(:,1)./x4(:,3);
x4(:,2) = x4(:,2)./x4(:,3);


% Assigning the output variables.
node4{1} = x1;
node4{2} = x2;
node4{3} = x3;
node4{4} = x4;

return

function [R] = tri10(xi,eta)
%----------------------------------------tri10---------------------------------%
% TRI10 is the finite element subroutine for a 10 node triangular element. It
% takes as inputs the control net of the triangle in physical space, B and the
% point in isoparametric space at which to evaluate the basis functions,
% [xi,eta]. It outputs the value of the basis function, R, its derivative with
% respect to physical coordinates, dR_dx, and the Jacobian determinate of the
% mapping from isoparametric space to physical space, J_det.


%------------------------------------------------------------------------------%
% Element parameters.
n = 3;
nen = 10;
% Find the barycentric coordinates of xi and eta
vert = [0,0;1,0;0,1];
x1 = vert(1,1);
x2 = vert(2,1);
x3 = vert(3,1);
y1 = vert(1,2);
y2 = vert(2,2);
y3 = vert(3,2);

detA =  det([x1,x2,x3; y1,y2,y3;1, 1, 1]);
detA1 = det([xi,x2,x3;eta,y2,y3;1, 1, 1]);
detA2 = det([x1,xi,x3;y1,eta,y3;1, 1, 1]);
detA3 = det([x1,x2,xi;y1,y2,eta;1, 1, 1]);

u = detA1/detA;
v = detA2/detA;
w = detA3/detA;

% Initializing variables
R = zeros(nen,1);

% Index and tuples. The index is the location in the control net (row,col) of
% the ith control point. Tuples is the index in barycentric coordinates of the
% ith control point.
idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];
tuples  = [ 3 0 0;...
    0 3 0;...
    0 0 3;...
    2 1 0;...
    1 2 0;...
    0 2 1;...
    0 1 2;...
    1 0 2;...
    2 0 1;...
    1 1 1];

% Loop through the control points.
for nn = 1:nen
    
    i = tuples(nn,1);
    j = tuples(nn,2);
    k = tuples(nn,3);
    
    % From page 141 of Bezier and B-splines. Calculate the ith basis function
    % its derivative with respect to barycentric coordinates.
    R(nn) = factorial(n)/...
        (factorial(i)*factorial(j)*factorial(k))*u^i*v^j*w^k;
    
    if i-1 < 0
        dR_du(nn,1) =0;
    else
        dR_du(nn,1) = n * factorial(n-1)/...
            (factorial(i-1)*factorial(j)*factorial(k))*u^(i-1)*v^j*w^k;
    end
    
    if j-1 < 0
        dR_du(nn,2) = 0;
    else
        dR_du(nn,2) = n * factorial(n-1)/...
            (factorial(i)*factorial(j-1)*factorial(k))*u^(i)*v^(j-1)*w^k;
    end
    
    if k-1 < 0;
        dR_du(nn,3) = 0;
    else
        dR_du(nn,3) = n * factorial(n-1)/...
            (factorial(i)*factorial(j)*factorial(k-1))*u^i*v^j*w^(k-1);
    end
    
end

return
