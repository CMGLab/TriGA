function [NODE,IEN, BFLAG,temp] = gambitFileIn(FileName)

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
for k = 1:K
    line   = fgetl(Fid);
    tmpcon = sscanf(line, '%lf');
    IEN(k,:) = tmpcon(4:13);
end

IEN = IEN';
% skip through material property section
for i=1:4
    line = fgetl(Fid);
end;

while isempty(intersect(line,'ENDOFSECTION'))
    line   = fgetl(Fid);
end

line = fgetl(Fid); line = fgetl(Fid);

% boundary codes (defined in Globals2D)
BFLAG = zeros(K,4);

% Read all the boundary conditions at the nodes
group =1;
while  line ~= -1 & ~strcmp(line(1:8),'TIMESTEP')
    bcprop = sscanf(line, '%s %u %u %u %u');
    bcprop = bcprop(end-3:end);
    NBELEM = bcprop(2);
    TYPE = bcprop(4);
    for bb = 1:NBELEM
        line = fgetl(Fid);
        belem = sscanf(line,'%u');
        BFLAG(belem(1),:) = [belem(1),belem(3),group,TYPE];
    end
    group = group+1;
    line = fgetl(Fid); line = fgetl(Fid); line = fgetl(Fid);
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