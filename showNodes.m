function [] = showNodes(node)
hold on
x = node(:,1);
y = node(:,2);

scatter(x,y);
a = [1:10]'; b = num2str(a); c = cellstr(b);
dx = 0.02; dy = 0.02; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);