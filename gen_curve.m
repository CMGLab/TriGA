function [C] = gen_curve(P,kv,n)
% Gen curve 
h = .1;
u  = 0:h:1;
k  = size(P,1);

C = zeros(length(u),2);

for j = 1:length(u)
    num = zeros(1,2);
den = 0;
    for i = 1:k
        num = num+ deBoor(i,n,u(j),kv)*P(i,3)*P(i,1:2);
        den = den+ deBoor(i,n,u(j),kv)*P(i,3);
    end
    C(j,:) = num/den;
end