%% Topoplot Padi
% clear all
% clear; close; clc;

%% signal
% y = 1;

%% etiquetas 10-20%
load labels
chan = labels;
load layout
[~,pos] = ismember(chan,labels);

%% CSD for the EEG layout
M1.lab = labels(pos);
M1.xy = data(pos,[4 5]);
MyTopo_fun(y,M1.xy,M1.lab)

figure
for k = 1:1
%     subplot(1,4,k)
    MyTopo_fun(y,M1.xy,M1.lab)
    axis square
end