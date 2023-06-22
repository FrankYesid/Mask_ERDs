%%Funcion para graficar los canales C3, Cz  y C4 de los ERDs, almacenados en estructura.
% Para la base de datos de BCICIV_2a
% Input:
%   r         ...  Contiene el ERD de ambas clases de cada uno de los canales.
%                  Se encuentra en la estructura ... r {fold}{class}.ERDS{ch}.erds
%                  Importante identificar que 
% erdcolormap .. Corresponde  
% Output:
%       ... Grafica correspondiente de las dos clases para los tres
%           canales centrales.
%          Variable que retorna la función de prepareData2.m del toolbox
%          Biosig3.6.
%   F. Y. Zapata C. 2019 
%%
function [maxi,mini]=Graficas_ERDs_channels2(r,erdcolormap,labels,labels2,name,sav,Dir_save)

% ------------- Filter bank -------------
f_low  = 4;    % frecuencia mínima.
f_high = 40;   % frecuencia máxima.
Window = 4;    % Tamaño de la banda de frecuencia.
Ovrlap = 2;    % Traslape de las bandas.
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';

mont = zeros(1,22);
mont([8,10,12]) = 1;
border = 0.1;  % Border around figure
border_plots = 0.01;  % Border around each plot
plot_area = 1 - 2 * border;
fold = 1;
clas = 1;
% Topographic layout
% if isfield(r{fold}{clas}, 'montage')

  montage = [0 0 0 0 1 33 34 0 0 0 0;...
                0 0 0 2 3 37 36 35 0 0 0;...
                0 7 6 5 4 38 39 40 41 42 0;...
                0 8 9 10 11 47 46 45 44 43 0;...
                0 15 14 13 12 48 49 50 51 52 0;...
                0 16 17 18 19 32 56 55 54 53 0;...
                24 23 22 21 20 31 57 58 59 60 61;...
                0 0 0 25 26 30 63 62 0 0 0;...
                0 0 0 0 27 29 64 0 0 0 0;...
                0 0 0 0 0 28 0 0 0 0 0];
    plot_index = find(montage' >= 1);
    n_rows = size(montage, 1);
    n_cols = size(montage, 2);
% else  % create default layout
%     plot_index = 1:length(r{fold}{clas}.ERDS);
%     n_cols = ceil(sqrt(length(r{fold}{clas}.ERDS)));
%     if (length(r{fold}{clas}.ERDS) > 2)
%         n_rows = n_cols;
%     else
%         n_rows = 1;
%     end;
% end;
% % n_cols = chan;
n_rows = 2;
i_width = plot_area / length(plot_index);% n_cols;  % Width of one subplot
i_height = plot_area / n_rows;  % Height of one subplot
font_size = 1/15;%1/32;  % Default normalized axes font size
f = figure;
set(f, 'PaperOrientation', 'landscape');
if ~exist('OCTAVE_VERSION','builtin')
    set(f, 'PaperType', 'A4');
end;
set(f, 'PaperUnits', 'centimeters');
set(f, 'PaperPosition', [1, 1, 27.7, 19]);
set(f, 'Color', [1 1 1]);
set(f, 'DefaultAxesFontUnits', 'normalized');
set(f, 'DefaultAxesFontSize', font_size);
plot_index =[8,10,12];
load('H:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
% Invert color map so that ERS is blue and ERD is red
colormap(erdcolormap);
a = cell(2, length(plot_index));  % Contains the axes of the ERDS subplots
for clas = 1:2
    i_width = plot_area / length(plot_index);% n_cols;  % Width of one subplot
    counter_total = 1;  %1 Iterates through all rows and columns
    if clas==1
        counter_plots = 1;  % Iterates through all subplots
        counter_plots11 = 1;
    else
        counter_plots =  length(plot_index)+1;  % Iterates through all subplots
        counter_plots11 = 1;
    end
    i_rows = clas;
    for i_cols = 1:length(plot_index)
        %         if sum(counter_total == plot_index) == 1
        %             i_rows = 1; i_cols = 1;
        a{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), (plot_area + border + border_plots - i_height * i_rows)+0.05, i_width - border_plots, i_height - border_plots]);
        set(f, 'CurrentAxes', a{counter_plots});
%         if strcmp(r{fold}{clas}.refmethod, 'absolute')
%             imagesc(r{fold}{clas}.t_plot, r{fold}{clas}.f_plot, r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds'); % r.ERDS{counter_plots}.erds'
%         else
%             [xx,yy]=find(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds== max(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(:)));
%             r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(xx,yy) = 1.5;
%             [xx,yy]=find(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds== min(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(:)));
%             r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(xx,yy) = -1;
            
            imagesc(0:1/250:7, 1:17, r{fold}{clas}{plot_index(counter_plots11)},[-1,1.5]);
%             maimagesc(0:1/250:7, 1:17, r{fold}{clas}{plot_index(counter_plots11)},[-1,1.5]);xi(fold,clas,counter_plots11) = max(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(:));
%             mini(fold,clas,counter_plots11)  = min(r{fold}{clas}.ERDS{plot_index(counter_plots11)}.erds(:));
%         end;
        %set(gca, 'Tag', num2str(counter_plots));
        set(a{counter_plots}, 'ydir', 'normal');
%         if counter_plots == 1
%             title('C$_{3}$','FontSize',15,'Interpreter','latex')
%         elseif counter_plots == 2
%             title('C$_{z}$','FontSize',15,'Interpreter','latex')
%         elseif counter_plots == 2
%             title('C$_{4}$','FontSize',15,'Interpreter','latex')
%         end
        if clas == 2 && labels2 == 0
%             if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
%                 set(a{counter_plots}, 'XTickLabel', '','FontSize',50,'Interpreter','latex');
%                 a{counter_plots}.FontSize = 50;
                a{counter_plots}.TickLabelInterpreter = 'latex';
                a{counter_plots}.XTick=([0:7]);
                a{counter_plots}.XTickLabel = num2str([0:7]');
%             end;
        else
            a{counter_plots}.XTickLabel = {};
        end
        
            if counter_plots == 1 || counter_plots == 4
                if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
                    set(a{counter_plots}, 'YTickLabel', '','FontSize',50,'Interpreter','latex');
                end;
                if labels ==1
                    a{counter_plots}.YTick=([1:2:17]);
                    a{counter_plots}.YTickLabel = num2str([filter_bank(a{counter_plots}.YTick,:)]);
                    a{counter_plots}.TickLength = [0.0 0.0001];
                else
                    a{counter_plots}.YTick=([1:2:17]);
                    a{counter_plots}.YTickLabel = '';
                end
            else
                a{counter_plots}.YTick={};
                a{counter_plots}.YTickLabel = {};
                a{counter_plots}.TickLength = [0.0 0.0001];
            end
%             if labels2 == 1
%                 a{counter_plots}.XTickLabel = '';
%             end
        
        % a{counter_plots}.FontSize = a{counter_plots}.FontSize+0.05;
        % Draw lines for reference interval and cue
        v = axis;
        %         line([r{fold}{clas}.ref(1), r{fold}{clas}.ref(1)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
        %         line([r{fold}{clas}.ref(2), r{fold}{clas}.ref(2)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
%         if isfield(r{fold}{clas}, 'cue')
            line([2, 2], [v(3), v(4)], 'Color', 'r','LineWidth',2);
%         end;
        a{counter_plots}.TickLabelInterpreter = 'latex';
        line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
        line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
        counter_plots = counter_plots + 1;
        counter_plots11 = counter_plots11+1;
        %         end;
        counter_total = counter_total + 1;
    end;
end;

% h=text(-16,33,'Frequency'); set(h,'Rotation',90);
% h1=text(-5,-1.2,'Time'); %set(h,'Rotation',90);
%                 axes('position', [border + border_plots,0.95 , 1 - 2 * (border + border_plots), border], 'visible', 'off');
%                 text(0.5, 0, r{fold}{clas}.heading, 'FontUnits', 'normalized', 'FontSize', 1/4, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'Interpreter', 'none', 'FontWeight', 'bold');
%  colorbar('Position',[0.905 0.16 0.02 0.79]);
if sav == 1
%     saveas(gca,['G:\Dropbox\ERD\Main_ERDS\Graficas\ERD_sub',num2str(name)],'png')
        saveas(gca,[Dir_save,'ERD_sub',num2str(name)],'epsc') 
          saveas(gca,[Dir_save,'ERD_sub',num2str(name)],'png') 
%         saveas(figs,) 
%         print('-depsc2', [Dir_save,'ERD_sub',num2str(name),'_class3_4.png')
end
