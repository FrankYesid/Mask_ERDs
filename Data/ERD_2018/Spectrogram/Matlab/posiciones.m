% para EEGLAB o topoplots
% load('posiciones.mat')
% para el topoplot
filename = fullfile('pos64channels.txt');
fileID = fopen(filename,'r');
% fileID = 6;
C = textscan(fileID,'%f %f32 %f32 %f32 %f32 %s');
M1.lab = C{6};
M1.xy = [C{2},C{3}];
load HeadModel.mat
% load('electrodesBCICIV1.mat')
% load('final.mat')
rel = 1:22; % para colocar mas importancia al electrodo.
sub = [64];
a=1;
s = 1;%[1,2,6,7]
%     subplot(1,4,a)
figure(s)
set(0,'DefaultFigureWindowStyle','docked')
sel = 1:sub(s); % que canales selecciono.
%     M1.xy = ;
%     M1.lab = electrodes;
drawnow
MyTopo_fun(rel,sel,M1.xy,M1.lab,[min(rel) max(rel)],0,0,200)
caxis([-1,1])
% colormap(gca,cmap);
title(['Subject ' num2str(s)])
axis square
axis off
a=a+1;


