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
function Graficas_ERDs_channels_all_giga(r,erdcolormap,name,sav,Dir_save,clas)
%% filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
%% Ubicación de los canales.
montage = '64ch';  % montaje de los canales que se van a utilizar.
channels = [ 1 33 34  2  3 37 36 35  7  6  5  4 38 39 40 41 42 ...
             8  9 10 11 47 46 45 44 43 15 14 13 12 48 49 50 51 ...
            52 16 17 18 19 32 56 55 54 53 24 23 22 21 20 31 57 ...
            58 59 60 61 25 26 30 63 62 27 29 64 28];   % número de canales que se van a utilizar.
% Toolbox relacionada
addpath(genpath('D:\Dropbox\ERD\results_ERDfc_subjects\pruebas_mask_ERDs\Toolbox\Biosig_ERD\biosig\t310_ERDSMaps\'))
%% Output parameters:
%   lap            ... Laplacian filter matrix.
%   plot_index ... Indices for plotting the montage.
%   n_rows     ... Number of rows of the montage.
%   n_cols      ... Number of columns of the montage.
fold = 1;                    % fold que se tiene en cuenta para la graficación.
[~, plot_index, n_rows, n_cols] = getMontage(montage);
border = 0.1;             % Border around figure
border_plots = 0.01;  % Border around each plot
plot_area = 1 - 2 * border;
i_width = plot_area / length(plot_index);% n_cols;  % Width of one subplot
i_height = plot_area / n_rows;              % Height of one subplot
font_size = 1/8;                                   %1/32;  % Default normalized axes font size
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
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
% Invert color map so that ERS is blue and ERD is red
colormap(erdcolormap);
i_width = plot_area/n_cols ; %length(plot_index);% n_cols;  % Width of one subplot
counter_total = 1;                   %1 Iterates through all rows and columns
counter_plots11 = 1;              %
a = cell(n_rows, n_cols);        % Contains the axes of the ERDS subplots
for i_rows = 1:n_rows
    counter_plots = 1;
    for i_cols = 1:n_cols
        if sum(counter_total == plot_index) == 1
            a{i_rows,i_cols} = axes('position', [border + border_plots + i_width * (i_cols - 1), (plot_area + border + border_plots - i_height * i_rows)+0.05, i_width - border_plots, i_height - border_plots]);
            set(f, 'CurrentAxes', a{i_rows,i_cols});
            %             if strcmp(r{fold}{clas}.refmethod, 'absolute')
            %                 imagesc(r{fold}{clas}.t_plot, r{fold}{clas}.f_plot, r{fold}{clas}.ERDS{counter_plots11}.erds'); % r.ERDS{counter_plots}.erds'
            %             else
            imagesc(-2:1/512:5,1:17, r{fold}{clas}{channels(counter_plots11)}, [-1, 1.5]);
            %             end;
            %set(gca, 'Tag', num2str(counter_plots));
            xlim([-1 4])
            set(a{i_rows,i_cols}, 'ydir', 'normal');
            if counter_plots == 1
                %             title('C$_{3}$','FontSize',15,'Interpreter','latex')
            elseif counter_plots == 2
                %             title('C$_{z}$','FontSize',15,'Interpreter','latex')
            elseif counter_plots == 3
                %             title('C$_{4}$','FontSize',15,'Interpreter','latex')
            end
            %             if r{fold}{clas}.class == 2
            %                 if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
            %             if
            %                     set(a{i_rows,i_cols}, 'XTickLabel', '','FontSize',50,'Interpreter','latex');
            %
            %                 end;
            %             else
            %                 a{i_rows,i_cols}.XTickLabel = {};
            %             end
            if channels(counter_plots11) == 14 || channels(counter_plots11) == 13 || channels(counter_plots11) == 12 || channels(counter_plots11) == 48 || channels(counter_plots11) == 49 || channels(counter_plots11) == 50 || channels(counter_plots11) == 51 %counter_total == 16 || counter_total == 18 || counter_total == 20
                a{i_rows,i_cols}.XColor = 'red';
                a{i_rows,i_cols}.YColor = 'red';
            end
            if  channels(counter_plots11) == 28 %counter_total == 4 || counter_total == 9 || counter_total == 15 || counter_total == 23 || counter_total == 31 || counter_total == 39
                if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
                    %                     set(a{i_rows,i_cols}, 'YTickLabel', '','FontSize',50,'Interpreter','latex');
                end;
%                                 if labels ==1
                a{i_rows,i_cols}.YTick=[1,3,5,7,9,11,13,15,17];
                a{i_rows,i_cols}.YTickLabel = num2str([filter_bank([1,3,5,7,9,11,13,15,17],:)]);
                a{i_rows,i_cols}.TickLength = [0.0 0.0001];
                %                 else
                %                     a{counter_plots}.YTick=r{fold}{clas}.f_plot([1,3,5,7,9,11,13]);
                %                     a{counter_plots}.YTickLabel = '';
                %                 end
            else
                a{i_rows,i_cols}.YTick={};
                a{i_rows,i_cols}.YTickLabel = {};
                a{i_rows,i_cols}.TickLength = [0.0 0.0001];
            end
            %             if labels2 == 1
            %                 a{counter_plots}.XTickLabel = '';
            %             end
            
            % a{counter_plots}.FontSize = a{counter_plots}.FontSize+0.05;
            % Draw lines for reference interval and cue
            v = axis;
            %         line([r{fold}{clas}.ref(1), r{fold}{clas}.ref(1)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
            %         line([r{fold}{clas}.ref(2), r{fold}{clas}.ref(2)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
%             if isfield(r{fold}{clas}, 'cue')
                line([0, 0], [v(3), v(4)], 'Color', 'r','LineWidth',2);
%             end;
            a{i_rows,i_cols}.TickLabelInterpreter = 'latex';
            line([0,0], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
            line([3,3], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
            if  channels(counter_plots11) ~= 28
                    a{i_rows,i_cols}.YTickLabel = {};
                   a{i_rows,i_cols}.XTickLabel = {};
            end
            counter_plots11 = counter_plots11+1;
        else
            %              a{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), (plot_area + border + border_plots - i_height * i_rows)+0.05, i_width - border_plots, i_height - border_plots]);
            %             set(f, 'CurrentAxes', a{counter_plots});
        end;
        counter_plots = counter_plots + 1;
        counter_total = counter_total + 1;
    end;
end;
% h=text(-16,33,'Frequency'); set(h,'Rotation',90);
% h1=text(-5,-1.2,'Time'); %set(h,'Rotation',90);
%                 axes('position', [border + border_plots,0.95 , 1 - 2 * (border + border_plots), border], 'visible', 'off');
%                 text(0.5, 0, r{fold}{clas}.heading, 'FontUnits', 'normalized', 'FontSize', 1/4, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'Interpreter', 'none', 'FontWeight', 'bold');
%  colorbar('Position',[0.905 0.16 0.02 0.79]);
if sav == 1
    %     saveas(gca,['G:\Dropbox\ERD\Main_ERDS\Graficas\Laplacian\ERD_sub',num2str(name)],'epsc')
    saveas(gca,[Dir_save,'ERD_sub',num2str(name),'class',num2str(clas),'_2'],'epsc')
    saveas(gca,[Dir_save,'ERD_sub',num2str(name),'class',num2str(clas),'_2'],'png')
end
