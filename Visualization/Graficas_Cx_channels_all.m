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
function [maxi,mini] = Graficas_ERDs_channels_all(r,erdcolormap,name,sav,Dir_save,clas)
%% Ubicación de los canales.
montage = '22ch'; % montaje de los canales que se van a utilizar.
channels = 1:22;   % número de canales que se van a utilizar.
% Toolbox relacionada al 
addpath(genpath('D:\Dropbox\ERD\Toolbox\Biosig_ERD\biosig\t310_ERDSMaps\'))
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
font_size = 1/15;                                   %1/32;  % Default normalized axes font size
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
% load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
% load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
% Invert color map so that ERS is blue and ERD is red
colormap(erdcolormap);
i_width = plot_area /  n_cols ; %length(plot_index);% n_cols;  % Width of one subplot
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
%             imagesc(r{fold}{clas}.t_plot, r{fold}{clas}.f_plot, r{fold}{clas}.ERDS{counter_plots11}.erds', [-1, 1.5]);
            imagesc(r{fold}{clas}.t_plot, r{fold}{clas}.f_plot, r{fold}{clas}.ERDS{counter_plots11}.erds',[0,1]);
            maxi(fold,clas,counter_plots11) = max(r{fold}{clas}.ERDS{counter_plots11}.erds(:));
            mini(fold,clas,counter_plots11)  = min(r{fold}{clas}.ERDS{counter_plots11}.erds(:));
            %             end;
            %set(gca, 'Tag', num2str(counter_plots));
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
            if counter_total == 16 || counter_total == 18 || counter_total == 20
                a{i_rows,i_cols}.XColor = 'red';
                a{i_rows,i_cols}.YColor = 'red';
            end
            if counter_total == 4 || counter_total == 9 || counter_total == 15 || counter_total == 23 || counter_total == 31 || counter_total == 39
                if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
                    %                     set(a{i_rows,i_cols}, 'YTickLabel', '','FontSize',50,'Interpreter','latex');
                end;
                %                 if labels ==1
                a{i_rows,i_cols}.YTick=r{fold}{clas}.f_plot([1,3,5,7,9,11,13,15,17]);
                a{i_rows,i_cols}.YTickLabel = num2str([r{1}{1}.f_low([1,3,5,7,9,11,13,15,17]);r{1}{1}.f_up([1,3,5,7,9,11,13,15,17])]');
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
            if isfield(r{fold}{clas}, 'cue')
                line([r{fold}{clas}.cue, r{fold}{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
            end;
            a{i_rows,i_cols}.TickLabelInterpreter = 'latex';
            line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
            line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
            if counter_total ~=  39
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
%     saveas(gca,[Dir_save,'ERD_sub',num2str(name),'class',num2str(clas),'_2'],'epsc')
%     saveas(gca,[Dir_save,'ERD_sub',num2str(name),'class',num2str(clas),'_2'],'png')
    matlab2tikz(['D:\Documento_tesis\Figures\Figs_cap3\Cx\maps\Map_Cx_all_bci2a_S_',num2str(name),'class',num2str(clas),'.tikz'])
end
