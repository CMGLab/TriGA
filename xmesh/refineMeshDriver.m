function [ NODE, IEN, BFLAG ] = refineMeshDriver( NODE, IEN, BFLAG, n )


% Convert to Edge based topological data structure.
TRI = reindex( IEN( 1:3, : )' );
TRI = TRI';
[EDGE, TRI2EDGE ] = IEN2EDGE(TRI);


for ii  = 1:n
    [ EDGE, TRI2EDGE, NODE, IEN, BFLAG ] = refineMeshTopo( EDGE, TRI2EDGE, NODE, IEN, BFLAG );
end

return

function [ EDGE, TRI2EDGE, NODE, IEN, BFLAG ] = refineMeshTopo( EDGEC, TRI2EDGEC, NODEC, IENC, BFLAGC )

side = [ 1 4 5 2;
        2 6 7 3; 
        3 8 9 1];
% Refine the topology
[EDGE, TRI2EDGE] = refineTopo( EDGEC, TRI2EDGEC );

nVert = max( EDGE( : ) );
nEdge = size( EDGE, 1 ); 
nelC = size( IENC, 2 );

% Initializing BFLAG variables. 
BFLAG = zeros(length(BFLAGC)*2,4);
ctr = 1;

IEN = zeros(10, 4*nelC);
NODE = zeros( nVert + 5*nEdge + 10* nelC, 3);
%figure, hold on
% Loop over the elements in the mesh.
for eeC = 1:nelC
    
    % Generate the local node array, and perform uniform refinement using
    % refine Triange.
    node = NODEC(IENC(:,eeC),:);
    node4 = refineTriangle(node);
    
    



    % Loop over sub elements and update connectivity.
    for rr = 1:4
    %    node4{rr};
    %    scatter(ans(:,1),ans(:,2))
        
        ee = 4*(eeC - 1) + rr;
        for ss = 1:3
            eidx = abs( TRI2EDGE( ee, ss ) );
            
            if TRI2EDGE( ee, ss ) > 0
                
                % Vertex
                IEN( side( ss, 1 ), ee ) = EDGE( eidx, 1 );
                NODE( EDGE( eidx, 1 ), : ) = node4{rr}( side( ss, 1 ), : );

                % Edge nodes.
                IEN( side( ss, 2 ), ee ) = nVert + 2*( eidx - 1 ) + 1 ;
                IEN( side( ss, 3 ), ee ) = nVert + 2*( eidx - 1 ) + 2 ;
                
            else
                % Vertex
                IEN( side( ss, 1 ), ee ) = EDGE( eidx, 2 );

                NODE( EDGE( eidx, 1 ), : ) = node4{rr}( side( ss, 4 ), : );
                
                % Edge nodes.
                IEN( side( ss, 2 ), ee ) = nVert + 2*( eidx - 1 ) + 2 ;
                IEN( side( ss, 3 ), ee ) = nVert + 2*( eidx - 1 ) + 1 ;
            end
            
            NODE( IEN( side( ss, 2:3 ), ee ), : ) = node4{rr}( side( ss, 2:3 ), : );
            
        end
        
        % Faces
        fidx = nVert + 2*nEdge + ee;
        IEN( 10, ee ) = fidx;
        NODE( fidx, : ) = node4{rr}( 10, : );
    end
    
    
    
    % Update the boundary condition matrix.
    temp  = BFLAG(BFLAG(:,1)==ee,:);
    for bb = 1:size(temp,1)
        for ss = 1:3
            if temp(bb,2) == ss
                BFLAG(2*(ctr-1)+1,:) = [4*(ee-1)+side( ss, 1 ) temp(bb,2:4)];
                BFLAG(2*(ctr-1)+2,:) = [4*(ee-1)+side( ss, 2 ) temp(bb,2:4)];
                ctr = ctr+1;
            end
        end
    end
    
end

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
    Rhat1(i,:) = tri10Bern(Xi1(i,1),Xi1(i,2))';
    Rhat2(i,:) = tri10Bern(Xi2(i,1),Xi2(i,2))';
    Rhat3(i,:) = tri10Bern(Xi3(i,1),Xi3(i,2))';
    Rhat4(i,:) = tri10Bern(Xi4(i,1),Xi4(i,2))';
    
end


for i = 1:10
    R(i,:) = tri10Bern(Xi(i,1),Xi(i,2))';
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

function [R] = tri10Bern(xi,eta)
%----------------------------------------tri10---------------------------------%
% TRI10BERN calculate the bernstein polynomials over the unit triangle.
%
% INPUT:
% xi: The xi location (in parametric space) at which to evaluate the basis
% functions.
%
% eta: The eta location (in parametric space) at which to evaluate the basis
% functions.
%
% OUTPUT:
% R: A 10x1 array containing the basis functions evaluated at [xi,eta].
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
end

return

function [EDGE4, TRI2EDGE4] = refineTopo(EDGE, TRI2EDGE)

nel = size(TRI2EDGE, 1);
nedge = size(EDGE, 1 );

TRI2EDGE4 = zeros(4*nel, 3);
EDGE4 = zeros( 2*nedge + 3*nel, 2);

ndof = max( EDGE(:) );

for ss = 1:nedge
    EDGE4( 2 * ( ss - 1 ) + 1, [ 1 2 ] ) = [ EDGE( ss, 1 ), (ndof + ss) ];
    EDGE4( 2 * ( ss - 1 ) + 2, [ 1 2 ] ) = [ (ndof + ss), EDGE( ss, 2 ) ];
end



% Loop over elements.
for ee = 1:nel
    
    % Loop over edges on the element.
            eloc  =zeros(9,1);
    for ss = 1:3
        eidx  = TRI2EDGE( ee, ss );
        
        if eidx > 0
            eloc( 2 * ( ss - 1) + 1 ) = 2 * ( abs( eidx) - 1) + 1 ;
            eloc( 2 * ( ss - 1) + 2 ) = 2 * ( abs( eidx) - 1) + 2 ;
        else
            eloc( 2 * ( ss - 1) + 1 ) = - ( 2 * ( abs( eidx) - 1) + 2 ) ;
            eloc( 2 * ( ss - 1) + 2 ) = - ( 2 * ( abs( eidx) - 1) + 1 ) ;
        end
    end
    eloc ( 7 ) = 2 * nedge + 3*(ee-1) + 1;
    eloc ( 8 ) = 2 * nedge + 3*(ee-1) + 2;
    eloc ( 9 ) = 2 * nedge + 3*(ee-1) + 3;
    
       
    
    EDGE4(eloc( 7 ) , [ 1 2] ) = [ ndof+abs(TRI2EDGE( ee, 2 )), ndof+abs(TRI2EDGE( ee, 3 )) ];
    EDGE4(eloc( 8 ), [ 1 2] ) = [ ndof+abs(TRI2EDGE( ee, 3 )), ndof+abs(TRI2EDGE( ee, 1 )) ];
    EDGE4(eloc( 9 ), [ 1 2] ) = [ ndof+abs(TRI2EDGE( ee, 1 )), ndof+abs(TRI2EDGE( ee, 2 )) ];



    TRI2EDGE4( 4*(ee-1)+ 1, [ 1 2 3 ] ) = [  eloc( 1 ), -eloc( 8 ),  eloc( 6 ) ];     
    TRI2EDGE4( 4*(ee-1)+ 2, [ 1 2 3 ] ) = [  eloc( 2 ),  eloc( 3 ), -eloc( 9 ) ];     
    TRI2EDGE4( 4*(ee-1)+ 3, [ 1 2 3 ] ) = [ -eloc( 7 ),  eloc( 4 ),  eloc( 5 ) ];     
    TRI2EDGE4( 4*(ee-1)+ 4, [ 1 2 3 ] ) = [  eloc( 7 ),  eloc( 8 ),  eloc( 9 ) ];     


end




return    
    
function TRI = reindex(TRI)

indices = unique( TRI(:) );

for ii = 1:length( indices ) 
   
    TRI(TRI == indices(ii) ) = ii;
    
    
end