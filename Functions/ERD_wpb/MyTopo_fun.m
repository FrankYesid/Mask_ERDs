function MyTopo_fun(Y,pos,label)

for i=1:2
    pos(:,i) = 0.9*((pos(:,i)-min(pos(:,i)))/range(pos(:,i))-0.5);
end
%figure
%HeadModel
load HeadModel
xc = HeadModel(1,:);
yc = HeadModel(2,:);
plot(xc,yc,'k')
hold on
%% Topoplot
GS = 100;
x = pos(:,1);
y = pos(:,2);
xi         = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
yi         = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
[Xi, Yi, Zi] = griddata(x', y, Y, xi', yi,'v4'); % interpolate the topographic data
%% Creating data mask
[TH R] = cart2pol(Xi,Yi);
Zi(R>0.5) = NaN;
deltax = xi(2)-xi(1); % length of grid entry
deltay = yi(2)-yi(1); % length of grid entry
h = surf(Xi-deltax/2, Yi-deltay/2, zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', 'flat');
hold on
contour(Xi-deltax/2, Yi-deltay/2,Zi,'k')
hold on
plot(x,y,'*')
set(gca,'XTick',[],'YTick',[]);
if nargin == 3
    text(x,y+0.02,label);
end
function [xunit yunit] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
