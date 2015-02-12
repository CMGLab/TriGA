function [p,t,kvloc,fnum,stats] = xmeshfaces(node,edge,P,KV,kvloc,face,hdata,options)
% ------------------------------------------------------------------------%
% XMESHFACES This function is a HEAVILY modified version of meshfaces, which
% is part of MESH2D (written by Darren Engwirda). 
%
% INPUT:
% node: an nx2 array of nodal coordinates that define the polygons that
% approximate the input NURBS curves. 
%
% edge: an ex2 array of edge connectivity information. 
%
% P: 1xnCurves cell array with each cell containing the control points of
% one of the input NURBS Curves. 
%
% KV: 1xnCurves cell array with each cell containing the knot vector of
% one of the input NURBS Curves. 
%
% kvloc: The knot locations corresponding to the points in node. 
%
% face: [] Right now this is a junk variable as xmesh only handles one
% face, but capabilities for more than one face should be added in the
% future. 
%
% hdata: [] Data struct containing user definied sizing data. Right now
% xmesh doesn't support this, so its also a junk variable. 
%
% options: Output options. If options.output is true, xmesh will give more
% information to the user. 
%
% OUTPUT: 
% p: an nx2 list of points giving the linear mesh that parameterizes the
% geonmetry. 
%
% t: an nelx3 array containing the triangulation connectivity information.
%
% kvloc: A refined kvloc vector that now corresponds to the refine polygons
% that make up the geometry boundary. 
%
% fnum: Number of faces (junk variable)
% stats: Mesh statistics
% ------------------------------------------------------------------------%

% Mesh timers. 
ts = cputime;
tic

% Get user options
options = getoptions(options);

% Dynamic quadtree decompositions. Call quadtree to generate the background
% mesh. Then discritise the boundary, and update the polygon that
% approximates our geometry. Keep doing this until boundarynodes does not
% split the input polygon, or until a maximum number of interations is
% reached.
nCurves = numel(P);
it = 1;
maxit = 10;
splitFlag = 1;

% Quadtree decomposition
%  PH    : Background mesh nodes
%  TH    : Background mesh triangles
%  HH    : Size function value at PH
[qtree.p,qtree.t,qtree.h] = quadtree(node,edge,hdata,options.dhmax,options.output);

while splitFlag && it <= maxit
    % Discretise edges
    [~,~,kvloc, splitFlag] = boundarynodes(qtree.p,qtree.t,qtree.h,node,edge,kvloc,options.output);
    
    % Find where the inserted polygonal node points lie om the actual curve.
    % Right now this methoda assumes a linear parameterization and tha each
    % side of the polygon is perfectly bisected. Will need to add
    % capability for boundary smoothing later.
    
    node = [];
    edge = [];
    ctr = 1;
    cStart = 1;
    for cc = 1:nCurves        
        kvtemp = unique( kvloc(kvloc>=cc-1 & kvloc<=cc) -(cc-1));
             
        for kk = 1:length(kvtemp)-1
            node(ctr,:)  = evalNURBS(P{cc},KV{cc},kvtemp(kk)); %#ok<*AGROW>
            edge(ctr,:) = [ctr ctr+1];
            ctr = ctr+1;
        end
        edge(ctr-1,2) = cStart;
        
        cStart = ctr;
    end
        
    it = it+1;
end

if it == maxit
    disp('Warning: meshfacescurved exited before converging')
    disp('Please email unmeshable geometries to engvall@colorado.edu')
end
pbnd = node;


% Mesh the connected face defined by node and edge. 
p = []; t = []; fnum = [];

face{1} = 1:length(pbnd);
for k = 1:length(face)
    
    % Mesh kth polygon
    [pnew,tnew] = meshpoly(node,edge(face{k},:),qtree,pbnd,options);
    
    % Add to global lists
    t = [t; tnew+size(p,1)];
    p = [p; pnew];
    fnum = [fnum; k*ones(size(tnew,1),1)];
    
end

% Ensure consistent, CCW ordered triangulation
[p,t,~,fnum] = fixmesh(p,t,[],fnum);

% Element quality
q = quality(p,t);

% Method statistics
stats = struct('Time',cputime-ts,'Triangles',size(t,1), ...
    'Nodes',size(p,1),'Mean_quality',mean(q),'Min_quality',min(q));

end      % meshfaces()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,e,kvloc, splitFlag]  = boundarynodes(ph,th,hh,node,edge,kvloc,output)

% Discretise the geometry based on the edge size requirements interpolated
% from the background mesh.

p = node;
e = edge;
i = mytsearch(ph(:,1),ph(:,2),th,p(:,1),p(:,2));
h = tinterp(ph,th,hh,p,i);

if output
    fprintf('Placing Boundary Nodes\n');
end

% Edge length
dxy = p(e(:,2),:)-p(e(:,1),:);
L = sqrt(sum(dxy.^2,2));
% Size function on edges
he = 0.5*(h(e(:,1))+h(e(:,2)));
% Split long edges
ratio = L./he;
split = (ratio>=1.5);

if any(split)
    splitFlag = 1;
    % Split edge at midpoint
    n1 = e(split,1);
    n2 = e(split,2);
    pm = 0.5*(p(n1,:)+p(n2,:));
    n3 = (1:size(pm,1))' + size(p,1);
    % New lists
    e(split,:) = [n1,n3];
    e = [e; n3,n2];
    p = [p; pm];
    
    % Size function at new nodes
    disp('tshearch')
    i = mytsearch(ph(:,1),ph(:,2),th,pm(:,1),pm(:,2));
    h = [h; tinterp(ph,th,hh,pm,i)];
    
    splitloc = find(split);
    for kk = 1:length(splitloc);
        addkv(kk) = (kvloc(splitloc(kk)) + kvloc(splitloc(kk)+1))/2;
    end
    kvloc = sort([kvloc,addkv]);
else
    splitFlag = 0;
    return
end
clc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = getoptions(options)

% Extract the user defined options

% Defaults
d_mlim   = 0.02;
d_maxit  = 20;
d_dhmax  = 0.3;
d_output = true;

if ~isempty(options)
    if ~isstruct(options)
        error('OPTIONS must be a structure array');
    end
    if numel(options)~=1
        error('Options cannot be an array of structures');
    end
    fields = fieldnames(options);
    names = {'mlim','maxit','dhmax','output'};
    for k = 1:length(fields)
        if strcmp(fields{k},names)
            error('Invalid field in OPTIONS');
        end
    end
    if isfield(options,'mlim')                                              % Movement tolerance
        checkposscalar(options.mlim,'options.mlim');
    else
        options.mlim = d_mlim;
    end
    if isfield(options,'maxit')                                             % Maximum iterations
        options.maxit = round(checkposscalar(options.maxit,'options.maxit'));
    else
        options.maxit = d_maxit;
    end
    if isfield(options,'dhmax')                                             % Size function gradient limit
        checkposscalar(options.dhmax,'options.dhmax');
    else
        options.dhmax = d_dhmax;
    end
    if isfield(options,'output')                                            % Output on/off
        checklogicalscalar(options.output,'options.output');
    else
        options.output = d_output;
    end
else                                                                       % Default values
    options.mlim   = d_mlim;
    options.maxit  = d_maxit;
    options.dhmax  = d_dhmax;
    options.output = d_output;
end
options.debug = false;

end      % getoptions()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checkposscalar(var,name)

% Helper function to check if var is a positive scalar.

if var<0 || any(size(var)>1)
    error([name,' must be a positive scalar']);
end

end      % checkposscalar()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checklogicalscalar(var,name)

% Helper function to check if var is a logical scalar.

if ~islogical(var) || any(size(var)>1)
    error([name,' must be a logical scalar']);
end

end      % checklogicalscalar()
