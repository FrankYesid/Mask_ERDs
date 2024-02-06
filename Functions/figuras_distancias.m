%%

% function Grafica_mask(mask)
clear; clc
% Definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';

% tamaño de la figura
border = 0.1;  % Border around figure
border_plots = 0.0015;  % Border around each plot
plot_area = 1 - 2 * border;
% Topographic layout
% if isfield(r{fold}{clas}, 'montage')
%     plot_index = find(r{fold}{clas}.montage' == 1);
%     n_rows = size(r{fold}{clas}.montage, 1);
%     n_cols = size(r{fold}{clas}.montage, 2);
% else  % create default layout
%     plot_index = 1:length(r{fold}{clas}.ERDS);
%     n_cols = ceil(sqrt(length(r{fold}{clas}.ERDS)));
%     if (length(r{fold}{clas}.ERDS) > 2)
%         n_rows = n_cols;
%     else
%         n_rows = 1;
%     end;
% end;
% n_cols = chan;

n_rows = 2;
plot_index = 1:5;
i_width = plot_area / length(plot_index);% n_cols;  % Width of one subplot
i_height = plot_area / n_rows;  % Height of one subplot
font_size = 1/32;  % Default normalized axes font size
for s = 1:9
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

% plot_index = find(mont==1);
% load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
% Invert color map so that ERS is blue and ERD is red
% colormap(erdcolormap);
% Carga el número de los sujetos.
% load('G:\Dropbox\ERD\Main_ERDS\Dis_corr500.mat')
 
fold = 1;
for clas = 1:2
    counter_total = 1;  %1 Iterates through all rows and columns
    %     if r{fold}{clas}.classes==1
    %         counter_plots = 1;  % Iterates through all subplots
    %         counter_plots11 = 1;
    %     else
        counter_plots = 1;  % Iterates through all subplots
        counter_plots11 = 1;
    %     end
    a = cell(2, length(plot_index));  % Contains the axes of the ERDS subplots
    i_rows = clas;
    for i_cols = 1:length(plot_index)
        %         if sum(counter_total == plot_index) == 1
        %             i_rows = 1; i_cols = 1;
        a{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), (plot_area + border + border_plots - i_height * i_rows)+0.05, i_width - border_plots, i_height - border_plots]);
        set(f, 'CurrentAxes', a{counter_plots});
        imagesc(mask_i{s,fold},[0 1]);
        %set(gca, 'Tag', num2str(counter_plots));
        set(a{counter_plots}, 'ydir', 'normal');
        %         if counter_plots == 1
        %             title('C$_{3}$','FontSize',8,'Interpreter','latex')
        %         elseif counter_plots == 2
        %             title('C$_{z}$','FontSize',8,'Interpreter','latex')
        %         elseif counter_plots == 3
        %             title('C$_{4}$','FontSize',8,'Interpreter','latex')
        %         end
       
        if counter_plots == 1 || counter_plots == 6
            %             if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
%             set(a{counter_plots}, 'YTickLabel', '','FontSize',30,'Interpreter','latex');
%             end;
%             a{counter_plots}.YTick=1;
            a{counter_plots}.YTickLabel = num2str(filter_bank);
            a{counter_plots}.TickLength = [0.0 0.0001];
            a{counter_plots}.TickLabelInterpreter = 'latex';
%             a{counter_plots}.FontSize = 30;
        else
            a{counter_plots}.YTick={};
            a{counter_plots}.YTickLabel = {};
            a{counter_plots}.TickLength = [0.0 0.0001];
            a{counter_plots}.TickLabelInterpreter = 'latex';
        end
         if clas == 2
            %             if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
            a{counter_plots}.XTick = 1:22;
%             set(a{counter_plots}, 'XTickLabel', '','FontSize',30,'Interpreter','latex');
            a{counter_plots}.XTickLabel = num2str([1:22]');
            a{counter_plots}.TickLength = [0.0 0.0001];
            a{counter_plots}.TickLabelInterpreter = 'latex';
            %                 end;
         else
             a{counter_plots}.XTickLabel = {};
         end
         % a{counter_plots}.FontSize = a{counter_plots}.FontSize+0.05;
         % Draw lines for reference interval and cue
         v = axis;
         %         line([r{fold}{clas}.ref(1), r{fold}{clas}.ref(1)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
         %         line([r{fold}{clas}.ref(2), r{fold}{clas}.ref(2)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
         %         if isfield(r{fold}{clas}, 'cue')
         %             line([r{fold}{clas}.cue, r{fold}{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
         %         end;
         %          a{counter_plots}.TickLabelInterpreter = 'latex';
         line([8, 8], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
         line([10, 10], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
         line([12, 12], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
         counter_plots = counter_plots + 1;
         fold =fold +1;
         counter_plots11 = counter_plots11+1;
         %         end;
        counter_total = counter_total + 1;
    end;
end;
h=text(-95,14,'Frequency','Interpreter','latex','FontSize',15); set(h,'Rotation',90);
h1=text(-36,-1.2,'Channels','Interpreter','latex','FontSize',15);
colorbar('Position',[0.905 0.15 0.02 0.8])
saveas(gcf,['G:\Dropbox\ERD\Main_ERDS\Graficas\Folds_organizada_500ms',filesep,'Subj_' num2str(s) '_Correlation '],'epsc')
% saveas(gcf,['G:\Dropbox\ERD\Main_ERDS\Graficas\Folds_organizada',filesep,'Subj_' num2str(s) '_Correlation '],'png')
matlab2tikz(['G:\Dropbox\ERD\Main_ERDS\Graficas\Folds_organizada_500ms',filesep,'Subj_' num2str(s) '_Correlation .tex'])
close
end;
