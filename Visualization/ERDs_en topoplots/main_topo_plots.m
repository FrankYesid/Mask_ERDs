% clear all; close all; clc
% Input:
%           rel         ... Es la matriz de 3 dimensiones (class x time x channels)
%           limits_x ... limites en el eje x.
%           limits_y ... limites en el eje y.
%           sel        ... canales seleccionados a graficar.
%           tam_label  ... tama�o de los labels.
%           wi         ... Ancho del plot.
%           hi          ... Alto del plot.
%           labels    ... (0 no)  o (1 si) quiere no los n�meros en los ejes.
% Output:
%           Figura...
% 
% 
function main_topo_plots(rel,limits_x,limits_y,sel,tam_label,wi,hi,labels)
% load positions electrodes
%% Filter bank
f_low  = 4; f_high = 40;
Window = 4; Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;f_low+Window:Ovrlap:f_high]'; %

%% time
time = 0:1/250:7;

%% positions
load('electrodesBCICIV2a.mat')
% load positions of plots
load('posiciones_plot.mat')
load('HeadModel.mat')    % model of the head.
Fil = 3; Col = 3;              % topoplot en subplots para cada uno de los sujetos.
t_e = 40;                        % tama�o visual de los electrodos seleccionados.
M1.xy = elec_pos(:,1:2);% posicion de los canales.
M1.lab = Channels;        % nombre de los canales.
pos = M1.xy;
label = M1.lab;
warning off
for i=1:2
    pos(:,i) = 0.9.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
end
xc = HeadModel(1,:);
yc = HeadModel(2,:);
hold on
% Topoplot
x = pos(:,1);
y = pos(:,2);
tmp = [x,y,x*0]*rotz(2);
x = tmp(:,1);
y = tmp(:,2);
tam = 1; % size of head
% x1 = 0.02;     y1 = 0.06;     w = 0.95;      h = 0.92;
% ax=axes('position',[x1,y1,w,h]);
plot(tam*xc,tam*yc,'color',[156,156,156]./255,'LineWidth',2)
axis off
hold on
% scatter(0.9.*x(sel),0.9.*y(sel),t_e,'b','filled')
load erdscolormap.mat
colormap(erdcolormap)
for ch = 1:numel(sel)
    ax1 = axes('Position',[posi(sel(ch),1) posi(sel(ch),2) wi hi]);
%     plot(ax1,squeeze(rel(:,:,sel(ch)))');

    imagesc(time,1:17,rel{ch},[-1,1.5])
    axis xy
%     xlim(limits_x); ylim(limits_y)
    %      set(ax1,'Gri
    set(ax1,'TickLabelInterpreter','latex','FontSize',tam_label,'XColor',[55,55,55]./255,'YColor',[55,55,55]./255,'YTick',1:3:17,'YTickLabel',num2str(filter_bank(1:3:17,:)))
    v = axis;
     line([0.5, 0.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
     line([2, 2], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
    line([2, 2], [v(3), v(4)],'LineStyle','--', 'Color', 'r');
    if labels == 0
        ax1.YTickLabel  = {}; ax1.XTickLabel  = {};
    end
    axis square
end
% set(gca,'axis','')
% 
% hold off
% axis square



