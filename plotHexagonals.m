function plotHexagonals(Y,ymax, xmax)
% Plots the solution Y on a hexagonal grid with dimensions
% YMAX by XMAX. 

Y = Y(:); 

for idx=1:length(Y)
  y = ymax - mod(idx-1,ymax);
  x = 1 + floor((idx-1)/(ymax)) - y/2;  
  coordinates = 2.*[x,y];
  cellPatch(coordinates,Y(idx));
end
axis off

function cellPatch(coordinates,fcolor)
% create hexagonal patch object
y = [0.5 1.5 0.5 -0.5 -1.5 -0.5 0.5];
x = [-1 0 1 1 0 -1 -1];
patch(x+coordinates(1),y+coordinates(2),fcolor, 'EdgeColor', 'black');
n = 50;

R = linspace(0.01,0.47,n);
B = linspace(0.01,0.19,n);
G = linspace(0.1,0.8,n);

colormap( [R(:), G(:), B(:)] );
