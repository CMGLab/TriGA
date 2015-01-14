function [tri] = cart2tri(n)
% cart2tri, calculates the vertices for a structured triangular grid over
% the triangular domain defined by and uper triangular nxn matrix.
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