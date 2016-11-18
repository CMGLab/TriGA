function [NODE,IEN, BFLAG,CFLAG,MAT,temp] = gambitFileInCAD(FileName)

Fid = fopen([FileName,'.neu'], 'rt');

% read intro
for i=1:6
    line = fgetl(Fid);
end

% Find number of nodes and number of elements
dims = fscanf(Fid, '%d');
Nv = dims(1); K = dims(2);

for i=1:2
    line = fgetl(Fid);
end

% read node coordinates
VX = (1:Nv); VY = (1:Nv); VW = (1:Nv);
for i = 1:Nv
    line = fgetl(Fid);
    tmpx = sscanf(line, '%lf');
    VX(i) = tmpx(2); VY(i) = tmpx(3); VW(i) = tmpx(4);
end

NODE = [VX',VY',VW'];
for i=1:2
    line = fgetl(Fid);
end

% read element to node connectivity
IEN = zeros(K, 10);
CFLAG = false(K,1);
for k = 1:K
    line   = fgetl(Fid);
    tmpcon = sscanf(line, '%lf');
    IEN(k,:) = tmpcon(4:13);
    if tmpcon(2)<0
        CFLAG(k) = true;
    end
end

IEN = IEN';

head = sscanf(line,'%s');
% Read material property section.

% Allocating a spot for the material array

MAT = zeros(0,2);
while line ~= -1 & ~strcmp(head(1:8),'BOUNDARY') %#ok<AND2>
    head = sscanf(line,'%s');

    % Read lines until we arrive at the next ELEMENT GROUP header.
    while line ~= -1 & ~strcmp(head(1:7),'ELEMENT')  %#ok<AND2>
        line = fgetl(Fid);
        head = sscanf(line,'%s');
    end
    
    % Read the next line, and get the group information.
    line = fgetl(Fid);
    material = sscanf(line,'%s %u %s %u %s %s');
    
    mat_grp = material(7);
    n_grp = material(16);
    
    matloc =  [ones(n_grp,1)*mat_grp zeros(n_grp,1)];
    
    line = fgetl(Fid); line = fgetl(Fid);
    
    n_line = floor(n_grp/10);
    n_last = mod(n_grp,10);
    
    ctr = 1;
    for ii = 1:n_line
        line = fgetl(Fid);
        mat_idx = sscanf(line,'%u');
        for jj = 1:10
            matloc(ctr,2) = mat_idx(jj);
            ctr = ctr+1;
        end
    end
    
    line = fgetl(Fid);
    mat_idx = sscanf(line,'%u');
    for jj = 1:n_last;
        matloc(ctr,2) = mat_idx(jj);
        ctr = ctr+1;
    end
    
    MAT = [MAT;matloc];
    
    line = fgetl(Fid);
    line = fgetl(Fid);
    head = sscanf(line,'%s');

    
end




BFLAG = zeros(K,4);

% Read all the boundary conditions at the nodes
group =1;
ctr = 1;
while  line ~= -1 & ~strcmp(line(1:8),'TIMESTEP')
    
    head = sscanf(line,'%s');

    % Read lines until we arrive at the next BOUNDARY GROUP header.
    while line ~= -1 & ~strcmp(head(1:8),'BOUNDARY')  %#ok<AND2>
        line = fgetl(Fid);
        head = sscanf(line,'%s');
    end
    
    % Read the next line that  has the group info 
    line = fgetl(Fid);
    bcprop = sscanf(line, '%s %u %u %u %u');
    bcprop = bcprop(end-3:end);
    NBELEM = bcprop(2);
    TYPE = bcprop(4);
    for bb = 1:NBELEM
        line = fgetl(Fid);
        belem = sscanf(line,'%u');
        BFLAG(ctr,:) = [belem(1),belem(3),group,TYPE];
        ctr = ctr+1;
    end
    group = group+1;
    line = fgetl(Fid);
    line = fgetl(Fid);
    if line ~= -1
        head = sscanf(line,'%s');
    end
end

line = fgetl(Fid);line = fgetl(Fid);

temp = zeros(size(NODE,1),1);
if line ~= -1
    for i = 1:Nv
        val = sscanf(line,'%f');
        temp(i)  = val(2);
        line = fgetl(Fid);
    end
end



% Close file
fclose(Fid);
return