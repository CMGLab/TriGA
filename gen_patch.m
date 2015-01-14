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

end

