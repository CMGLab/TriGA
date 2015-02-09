function [] = showL2(filename,elemL2)
% -----------------------------------------------------------------------------%
% SHOWRESULTS reads in the data from fea2d stored in 'fileName.mat' and prints a
% heatmap of the results.
% -----------------------------------------------------------------------------%
[NODE,IEN,~,temp] = gambitFileIn(filename);
nel = size(IEN,2);
nen = size(IEN,1);

idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];
for ee  = 1:nel
    node = NODE(IEN(:,ee),:);
    for i = 1:10
        B{ee}(idx(i,1),idx(i,2),1) = node(i,1);
        B{ee}(idx(i,1),idx(i,2),2) = node(i,2);
    end
    
end


figure
hold on
    for ee = 1:nel
        if mod(ee,100) == 0
        clc
            fprintf('Plotting element %4.0f of %4.0f \n',[ee,nel])
        end
        gen_patch(B{ee},elemL2(ee)*ones(10,1),1/3)
    end
    
axis equal 
colorbar
    
    return
    
        
function [] = gen_patch( B ,data,h)
% ---------------------------------------------------------------------------- %



idx = [1 1; 4 1; 1 4; 2 1; 3 1; 3 2; 2 3; 1 3; 1 2; 2 2];
for i  = 1:10
    B(idx(i,1),idx(i,2),4) = data(i);
end

u = 0:h:1;
v = 0:h:1;
n = 3;
% Preallocate the triangular surface coordinate matrices.
t1  =zeros(length(u),length(v));
t2  =zeros(length(u),length(v));
t3  =zeros(length(u),length(v));

% Calculate the triangular surface coordinate matrices.
for a = 1:length(u)
    for b =1:length(v)
        if u(a)+v(b)<=1
            for i = 0:n
                for j = 0:n
                    if (i+j)<=n
                        
                        binco = factorial(n)./(factorial(i)*factorial(j)*...
                            factorial(n-i-j));

                        coeff = binco*u(a)^i*v(b)^j*(1-u(a)-v(b))^(n-i-j);
                        t1(a,b) = t1(a,b)+ B(i+1,j+1,1)*coeff;
                        t2(a,b) = t2(a,b)+ B(i+1,j+1,2)*coeff;
                        t3(a,b) = t3(a,b)+ B(i+1,j+1,4)*coeff;
                    end
                end
            end
        end
    end
end
tri = cart2tri(length(u));
trisurf(tri,t1,t2,t3,'EdgeAlpha',0)

return

function [tri] = cart2tri(n)
% cart2tri, calculates the vertices for a structured triangular grid over
% the triangular domain defined by an upper triangular nxn matrix.
% Counterclockwise numbering starting with the top left vertex.

count = 1;
trirows = sum(1:n-1)+sum(1:n-2);
tri = zeros(trirows,3);
for c = 1:n-1
    for r = 1:n-c
        tri(count,1) = n*(c-1)+r;
        tri(count,2) = n*(c-1)+r+1;
        tri(count,3) = n*c+r;
        count = count+1;
        
        if r<n-c
            tri(count,1) = n*(c-1)+r+1;
            tri(count,2) = n*c+r+1;
            tri(count,3) = n*c+r;
            count = count+1;
        end
        
    end
end
    
