fr1 = 3;
fr2 = 3;
for ch = 1:size(Dis,1)    % tamaño de los canales.
    aa = 1;
    for fre = [3]               % tamaño de los
        for win = 1:size(Dis,3)
            if fr1 == fr2
                Dis2(ch,win) = squeeze(Dis(ch,fre,win));
            else
                Dis2(ch,aa,win) = Dis(ch,fre,win)./sum(Dis(ch,fr1:fr2,win));
            end
        end
        aa = aa+1;
    end
end
if fr1 == fr2
    Dis_ = Dis2;
else
    Dis_ = squeeze(sum(Dis2,2));
end





%         end
%         save(['G:\Dropbox\ERD\ERD_2019\BCI\BCICIV_2a_0' num2str(ss) filesep  num2str(ss) '_' num2str(clas) '_' num2str(rep) 'ERD.mat'],'erd')
%         save(['F:\BCI\BCICIV_2a_0' num2str(ss) filesep  num2str(ss) '_' num2str(clas) '_'  'ERD_g.mat'],'erd')

%% Graficar los datos relacionados al ERD
% % r = erd;
% % load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
% % % Invert color map so that ERS is blue and ERD is red
% % figure;
% % clas = 1;
% % subplot(2,2,1) 
% % fig = imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{8}.erds',[-1,1.5]); 
% % colormap(erdcolormap); axis xy; fig.Parent.TickLabelInterpreter='latex';
% % v = axis; line([r{clas}.cue, r{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
% % line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % 
% % subplot(2,2,2) 
% % imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{12}.erds',[-1,1.5]); 
% % colormap(erdcolormap); axis xy; fig.Parent.TickLabelInterpreter='latex';
% % v = axis; line([r{clas}.cue, r{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
% % line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % 
% % clas = 2;
% % subplot(2,2,3) 
% % imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{8}.erds',[-1,1.5]); 
% % colormap(erdcolormap); axis xy; fig.Parent.TickLabelInterpreter='latex';
% % v = axis; line([r{clas}.cue, r{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
% % line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % 
% % subplot(2,2,4) 
% % imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{12}.erds',[-1,1.5]); 
% % colormap(erdcolormap); axis xy; fig.Parent.TickLabelInterpreter='latex';
% % v = axis; line([r{clas}.cue, r{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
% % line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% % 
% % %%
% % border = 0.1;  % Border around figure
% % border_plots = 0.01;  % Border around each plot
% % plot_area = 1 - 2 * border;
% % % Topographic layout 
% % if isfield(r{clas}, 'montage')
% %     plot_index = find(r{clas}.montage' == 1);
% %     n_rows = size(r{clas}.montage, 1);
% %     n_cols = size(r{clas}.montage, 2);
% % else  % create default layout
% %     plot_index = 1:length(r{clas}.ERDS);
% %     n_cols = ceil(sqrt(length(r{clas}.ERDS)));
% %     if (length(r{clas}.ERDS) > 2)
% %         n_rows = n_cols;
% %     else
% %         n_rows = 1;
% %     end;
% % end;
% % % n_cols = chan;
% % n_rows = 2;
% % i_width = plot_area / length(plot_index);% n_cols;  % Width of one subplot
% % i_height = plot_area / n_rows;  % Height of one subplot
% % font_size = 1/32;  % Default normalized axes font size
% % f = figure;
% % set(f, 'PaperOrientation', 'landscape');
% % if ~exist('OCTAVE_VERSION','builtin')
% %     set(f, 'PaperType', 'A4');
% % end;
% % set(f, 'PaperUnits', 'centimeters');
% % set(f, 'PaperPosition', [1, 1, 27.7, 19]);
% % set(f, 'Color', [1 1 1]);
% % set(f, 'DefaultAxesFontUnits', 'normalized');
% % set(f, 'DefaultAxesFontSize', font_size);
% % plot_index = find(mont==1);
% % 
% % 
% % %%
% % for clas = 1:2
% %     counter_total = 1;  %1 Iterates through all rows and columns
% %     if r{clas}.classes==1
% %         counter_plots = 1;  % Iterates through all subplots
% %         counter_plots11 = 1;
% %     else
% %         counter_plots =  length(plot_index)+1;  % Iterates through all subplots
% %         counter_plots11 = 1;
% %     end
% %     a = cell(2, length(plot_index));  % Contains the axes of the ERDS subplots
% %     i_rows = clas;
% %     for i_cols = 1:length(plot_index)
% %         %         if sum(counter_total == plot_index) == 1
% %         %             i_rows = 1; i_cols = 1;
% %         a{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), (plot_area + border + border_plots - i_height * i_rows)+0.05, i_width - border_plots, i_height - border_plots]);
% %         set(f, 'CurrentAxes', a{counter_plots});
% %         if strcmp(r{clas}.refmethod, 'absolute')
% %             imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{plot_index(counter_plots11)}.erds'); % r.ERDS{counter_plots}.erds'
% %         else
% %             imagesc(r{clas}.t_plot, r{clas}.f_plot, r{clas}.ERDS{plot_index(counter_plots11)}.erds', [-1, 1.5]);
% %         end;
% %         %set(gca, 'Tag', num2str(counter_plots));
% %         set(a{counter_plots}, 'ydir', 'normal');
% %         if counter_plots == 1
% %             title('C$_{3}$','FontSize',8,'Interpreter','latex')
% %         elseif counter_plots == 2
% %             title('C$_{z}$','FontSize',8,'Interpreter','latex')
% %         elseif counter_plots == 3
% %             title('C$_{4}$','FontSize',8,'Interpreter','latex')
% %         end
% %         if r{clas}.class == 2
% %             if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
% %                 set(a{counter_plots}, 'XTickLabel', '','FontSize',30,'Interpreter','latex');
% %             end;
% %         else
% %             a{counter_plots}.XTickLabel = {};
% %         end
% %         if counter_plots == 1 || counter_plots == 4
% %             if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
% %                 set(a{counter_plots}, 'YTickLabel', '','FontSize',30,'Interpreter','latex');
% %             end;
% %             a{counter_plots}.YTick=r{clas}.f_plot;
% %             a{counter_plots}.YTickLabel = r{clas}.f_plot;
% %             a{counter_plots}.TickLength = [0.0 0.0001];
% %         else
% %             a{counter_plots}.YTick={};
% %             a{counter_plots}.YTickLabel = {};
% %             a{counter_plots}.TickLength = [0.0 0.0001];
% %         end
% %         % a{counter_plots}.FontSize = a{counter_plots}.FontSize+0.05;
% %         % Draw lines for reference interval and cue
% %         v = axis;
% %         line([r{clas}.ref(1), r{clas}.ref(1)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
% %         line([r{clas}.ref(2), r{clas}.ref(2)], [v(3), v(4)], 'LineStyle', '-.', 'Color', 'g','LineWidth',2);
% %         if isfield(r{clas}, 'cue')
% %             line([r{clas}.cue, r{clas}.cue], [v(3), v(4)], 'Color', 'r','LineWidth',2);
% %         end;
% %         a{counter_plots}.TickLabelInterpreter = 'latex';
% %         line([2.5, 2.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% %         line([4.5, 4.5], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k','LineWidth',1);
% %         counter_plots = counter_plots + 1;
% %         counter_plots11 = counter_plots11+1;
% %         %         end;
% %         counter_total = counter_total + 1;
% %     end;
% % end;
% % h=text(-16,33,'Frequency'); set(h,'Rotation',90);
% % h1=text(-5,-1.2,'Time'); %set(h,'Rotation',90);

%%
%     for sub = 1:9
%         set(0,'DefaultFigureWindowStyle','docked')
%         %         clear ERD ERD_ Dis_Eu ERD_plot ERD_mean Dis_Eu Dis_Eu_
%         load(['F:\New_ERD_abri_2019\ERD_folds_s' num2str(sub) '.mat'])
%         %         CL = [1,2];
%         %         for folds = 1:Nfolds
%         %             for cl = CL
%         %                 for ch = 1:22
%         %                     ERD{cl}(ch,:,:)= erd{folds}{cl}.ERDS{ch}.erds;
%         %                 end
%         %             end
%         % %         end
%
%         %         load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
%
%         %         for cl = CL
%         %             figure
%         %             imagesc(r1{1}{1}.t_plot,1:22,squeeze(ERD_{cl}(:,:,3)),[-1 1.5])
%         %             axis xy
%         %             colormap(erdcolormap);
%         %             colorbar
%         %             title(['Sub: ' num2str(sub) ' L: ' num2str(cl)])
%         %         end
%
%         %         for cl = CL
%         %             for ch = [8,10,12]
%         %                 figure
%         %                 imagesc(r1{1}{1}.t_plot,1:17,squeeze(ERD_{cl}(ch,:,:))',[-1 1.5])
%         %                 axis xy; fig = gca;
%         %
%         %                 colormap(erdcolormap);
%         %                 colorbar
%         %                 title(['Sub: ' num2str(sub) ' L:' num2str(cl) ' Ch:' num2str(ch)])
%         %             end
%         %         end
%         %         cl =2;
%         %         load(['S:\Mi unidad\Conecctividad funcional\window\WPLI_cx_new_2_all' num2str(sub)  num2str(cl) '_0_.mat'])
%         %         Cx_ = squeeze(Con_wpli_cwt);
%         %         if sub==8
%         %             % sub = 8
%         %             Cl_s_ = [[7,8];[7,12];[7,13];[7,14];[8,12];[8,13];[8,14];[12,13];[12,14];[13,14]];
%         %         else
%         %             % sub = 2
%         %             Cl_s_ = [[7,8];[7,12];[7,13];[7,14];[8,12];[8,13];[8,14];[12,13];[12,14];[13,14]];
%         %         end
%         %         load(['S:\Mi unidad\Conecctividad funcional\window\WPLI_frs_.mat'])
%         %         for chan = 1:size(Cl_s_,1)
%         % %             for cl = 1
%         %                 figure
%         %                 Cx = Cl_s_(chan,:);
%         %                 contourf(squeeze(Cx_(Cx(2),Cx(1),:,:)),'LineStyle','none')
%         %                 caxis([0 1]); fig = gca;
%         %                 colorbar
%         %                 set(gca,'YScale','log')
%         %                 fig.YTick = [1,10,20,30,40,50,60,70];
%         %                 fig.XTickLabel = {num2str(r1{1}{1}.t_plot(fig.XTick)','%2.1f')};
%         %                 fig.YTickLabel = {num2str(frs(fig.YTick)')};
%         %
%         %                 title(['Sub: ' num2str(sub) ' L:' num2str(cl) ' Chs:' num2str(Cx(1)) '-' num2str(Cx(2))])
%         % %             end
%         %         end
%         for folds = 1:Nfolds
%             for clase = 1:2
%                 for ch = 1:size(erd{folds}{clase}.ERDS,2)
%                     ERD{folds}{clase}{ch} = erd{folds}{clase}.ERDS{ch}.erds;
%                 end
%             end
%         end
%         %
%         %         % Ventanas de los ERDs.
%         %         bb = 1;
%         for len = 0.5
%             over = 0.1;
%             fs = r1{1}{1}.fs;
%             time = numel(r1{1}{1}.t_plot);
%             t = 1:round(len*fs*over):(time-(len*fs));
%             Npar1 = numel(t);
%             for folds = 1:Nfolds
%                 for cl = 1:2
%                     for ch = 1:size(erd{folds}{cl}.ERDS,2)
%                         for fre = 1:numel(r1{1}{1}.f_plot)
%                             for tao = 1:Npar1
%                                 ERD_{folds}{cl}{ch}{fre}{tao} = ERD{folds}{cl}{ch}(t(tao):(t(tao)+len*fs-1),fre);
%                             end
%                         end
%                     end
%                 end
%             end
%             t_ = t./250;
%             w = find((t_>2.5 & t_<4.5)==1);
%             Npar = numel(w);
%
%             % distancia euclidea promediando los N. folds
%             for folds = 1:Nfolds
%                 for ch = 1:size(erd{folds}{cl}.ERDS,2)
%                     for fre = 1:numel(r1{1}{1}.f_plot)
%                         for tao = 1:Npar
%                             Dis_Eu{folds}(ch,fre,tao) = ERD_{folds}{1}{ch}{fre}{w(tao)}'*ERD_{folds}{2}{ch}{fre}{w(tao)}; %norm(ERD_{folds}{1}{ch}{fre}{tao}-ERD_{folds}{2}{ch}{fre}{tao});
%                         end
%                     end
%                 end
%             end
%
%             %             for fold =1:Nfolds
%             %                 figure
%             %                 plot(squeeze(Dis_Eu{fold}(1:22,3,:))')
%             %                 fig = gca;
%             %                 fig.XTick = [1 numel(w)];
%             %                 xlim([1 numel(w)])
%             %                 fig.XTickLabel = {'2.5','4.5'};
%             %                 title(['Sub: ' num2str(sub) ' F: ' num2str(fold)])
%             %             end
%
%             % promedio de las ventanas
%             for fold = 1:Nfolds
%                 for ch = 1:size(erd{folds}{cl}.ERDS,2)
%                     for fre = 1:numel(r1{1}{1}.f_plot)
%                         ERD_mean_w{fold} = squeeze(mean(Dis_Eu{fold},3));
%                     end
%                 end
%             end
%
%
%             % plot de las ventanas de la correlacion del ERD en los folds
%             %             for fold =1:Nfolds
%             %                 figure
%             %                 imagesc(1:22,1:17,abs(ERD_mean_w{fold}(1:22,:))'./max(max(ERD_mean_w{fold}(1:22,:))))
%             %                 axis xy; fig = gca;
%             %                 fig.YTickLabel = num2str([r1{1}{1}.f_low(fig.YTick)', r1{1}{1}.f_up(fig.YTick)']);
%             %                 colorbar
%             %                 title(['Sub: ' num2str(sub) ' Fold:' num2str(fold) ' w: ' num2str(len) ])
%             %                 fig = gca;
%             %                 saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\dis_ERD_corr' filesep 'Subj_' num2str(sub) 'fold' num2str(fold)],'png')
%             %                 close
%             %                 %                 saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\dis_ERD_corr' filesep 'Subj_' num2str(sub) 'fold' num2str(fold)],'epsc')
%             %             end
%             for fold =1:Nfolds
%                 col = [215,215,215]; col2 = [014,041,075]; %oscuro
%                 c = [linspace(col2(1),col(1),64)',linspace(col2(2),col(2),64)',...
%                     linspace(col2(3),col(3),64)']/255;
%
%                 %             if sum(isnan(rho))>=1
%                 %                 rho(isnan(rho)==1) = threshold(end);
%                 %             end %para los nan detectados
%                 mask_i = abs(ERD_mean_w{fold}(1:22,:))'./max(max(ERD_mean_w{fold}(1:22,:)));
%                 s = sub;
%                 % promedio correlaciones by channel
%                 Pi_c = zeros(1,22);
%                 for i=1:22; Pi_c(i) = mean(mask_i(:,i));end
%
%                 % promedio correlaciones by frequency
%                 Pi_f = zeros(1,size(filter_bank,1));
%                 for i=1:17; Pi_f(i) = mean(mask_i(i,:)); end
%                 graf=1;
%                 %% grafica
%                 if graf == 1
%                     figure; x = 0.25;     y1 = 0.15;     w = 0.55;      h = 0.8;
%                     subplot('position',[x,y1,w,h]); imagesc(mask_i,[0,1]); %min(rho),max(rho)
%                     title(num2str(s),'FontSize',15,'Interpreter','latex')
%                     %         colormap(c)
%                     axis xy;  xticks([]);
%                     xcol = x+w+0.003+0.01+0.06+0.01;
%                     colorbar('Position',[xcol y1 0.01 h]);
%
%                     yticks(1:size(filter_bank,1));
%                     yticklabels({num2str(filter_bank)})
%                     ylabel('Frequency','FontSize',36,'Interpreter',  'latex')
%
%                     %     x = 0.25; y2 = 0.12;  w = 0.55; h1 = 0.10;
%                     %     subplot('position',[x,y2,w,h1])
%                     %     bar(1:22,Pi_c,'FaceColor',c(end,:),...
%                     %         'EdgeColor',c(end,:)); xlim([0.6,22.4]);ylim([0,1])
%                     xticks(1:22);
%                     xticklabels({1:7,'C3',9,'Cz',11,'C4',13:22})
%                     xlabel('Channels','FontSize',36,'Interpreter','latex')
%
%                     % marginal en el eje y
%                     x_m2 = x+w+0.003+0.01; y_3 = y1;%%%%
%                     w_m2 = 0.06; h_m2 = h;
%                     subplot('position',[x_m2,y_3,w_m2,h_m2])
%                     barh(1:17,Pi_f,'FaceColor',c(end,:),'EdgeColor',c(end,:));
%                     ylim([0.6,17.4]);xlim([0,1])
%                     yticks([])
%                     %     saveas(gcf,[SUBJECTS_DIR2 filesep 'masks' filesep SUBJECTS{s}],'epsc')
%                     %     saveas(gcf,[SUBJECTS_DIR2 filesep 'masks' filesep SUBJECTS{s}],'fig')
%                 end
%                 fig = gca;
%                 saveas(fig,['G:\Dropbox\ERD\results_ERDfc_subjects\dis_ERD_corr' filesep 'Subj_' num2str(sub) 'fold' num2str(fold)],'png')
%                 close
%             end
%
%
%
%             %             close all
%             %             % promedio de los folds.
%             %             for ch = 1:size(erd{folds}{cl}.ERDS,2)
%             %                 Dis_ERD_mean{ch} = squeeze(mean(Dis_Eu{ch},1));
%             %             end
%             %             % Promedio de los folds para plotear.
%             %             for folds = 1:Nfolds
%             %                 for clase = 1:2
%             %                     for ch = 1:size(erd{folds}{clase}.ERDS,2)
%             %                         ERD_plot{clase}{ch}(folds,:,:) = erd{folds}{clase}.ERDS{ch}.erds;
%             %                     end
%             %                 end
%             %             end
%             %
%             %             for clase = 1:2
%             %                 for ch =1:size(erd{folds}{clase}.ERDS,2)
%             %                     ERD_mean{clase}{ch} = squeeze(mean(ERD_plot{clase}{ch},1));
%             %                 end
%             %             end
%             %             set(0,'DefaultFigureWindowStyle','docked')
%             %             %Graficar las imagenes de ERD
%             %             a = 1;
%             %             load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
%             %
%             %             if len ==2
%             %                 figure(1)
%             %                 for clase = 1:2
%             %                     for ch = [8,10,12] %1:size(erd{folds}{clase}.ERDS,2)
%             %                         subplot(2,3,a)
%             %                         %             figure
%             %                         imagesc(r1{1}{1}.t_plot,r1{1}{1}.f_plot,ERD_mean{clase}{ch}',[-1 1.5])
%             %                         v = axis; fig = gca;
%             %                         line([r1{1}{1}.ref(1), r1{1}{1}.ref(1)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
%             %                         line([r1{1}{1}.ref(2), r1{1}{1}.ref(2)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
%             %                         if isfield(r1{1}{1}, 'cue')
%             %                             line([r1{1}{1}.cue, r1{1}{1}.cue], [v(3), v(4)], 'Color', 'k');
%             %                         end;
%             %                         fig.TickLabelInterpreter = 'latex';
%             %                         fig.YTick=r1{1}{1}.f_plot;
%             %                         fig.YTickLabel = r1{1}{1}.f_plot;
%             %                         axis xy; colorbar;
%             %                         colormap(erdcolormap);
%             %                         if ch == 8 && clase == 1
%             %                            title('C3')
%             %                         elseif ch == 10 && clase == 1
%             %                             title('Cz')
%             %                         elseif ch == 12 && clase == 1
%             %                             title('C4')
%             %                         end
%             %                         a = a+1;
%             %                     end
%             %                 end
%             %                 fig = gca;
%             %                 saveas(fig,['sub_' num2str(sub) 'erd.png'])%    saveas(fig,[''],'epsc')
%             %                 saveas(fig,['sub_' num2str(sub) 'erd.fig'])
%             %                 close
%             %                 %Distancia Euclidea
%             %             end
%             %
%             %             for ch =1:size(erd{folds}{clase}.ERDS,2)
%             %                 for freq = 1:size(ERD_mean{1}{ch},2)
%             %                     Dis_Eu{ch}(:,freq) = norm(ERD_mean{1}{ch}(:,freq) - ERD_mean{2}{ch}(:,freq));
%             %                 end
%             %             end
%             %
%             %             figure(2)
%             %             for ch = [8,10,12] %1:size(erd{folds}{clase}.ERDS,2)
%             %                 subplot(4,3,bb)
%             %                 %         figure
%             %                 for c = 1:22
%             %                     abc(c,:,:)=Dis_ERD_mean{c};
%             %                 end
%             %                 [~,b] = min(abs(t/250-6));
%             %                 imagesc(t./fs,1:17,Dis_ERD_mean{ch}./max(max(max(abc(:,:,1:b)))),[0 1])
%             %                 v = axis; fig = gca;
%             %                 line([r1{1}{1}.ref(1), r1{1}{1}.ref(1)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
%             %                 line([r1{1}{1}.ref(2), r1{1}{1}.ref(2)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
%             %                 if isfield(r1{1}{1}, 'cue')
%             %                     line([r1{1}{1}.cue, r1{1}{1}.cue], [v(3), v(4)], 'Color', 'k');
%             %                 end;
%             %                 if ch == 8 && len == 2
%             %                     title('C3')
%             %                     ylabel(['w ',num2str(len)])
%             %                 elseif ch == 10 && len == 2
%             %                     title('Cz')
%             %                 elseif ch == 12 && len == 2
%             %                     title('C4')
%             %                 end
%             %                 fig.TickLabelInterpreter = 'latex';
%             %                 fig.YTick = 1:17;
%             %                 fig.YTickLabel = r1{1}{1}.f_plot;
%             %                 axis xy; xlim([0 7]); colorbar
%             %
%             %                 %         colormap(erdcolormap);
%             %                 bb = bb+1;
%             %             end
%             %             clear abc
%             %         end
%             %         fig = gca;
%             %         saveas(fig,['sub_' num2str(sub) 'dis'],'fig')%    saveas(fig,[''],'epsc')
%             %         saveas(fig,['sub_' num2str(sub) 'dis'],'png')
%             %         close
%         end
%     end

%%    imagen de los canales en las diferentes frecuencias.
%     for ch = 1:size(erd{folds}{clase}.ERDS,2)
%         Dis_Eu_cf(ch,:,:) = Dis_ERD_mean{ch};
%     end
%     figure
%     imagesc(t./fs,1:17,Dis_Eu_cf(:,:,14)./max(max(cell2mat(Dis_ERD_mean([1:22])))),[0 1])
%     axis xy
%
%
%     % topoplot
%     for ch = 1:size(erd{folds}{clase}.ERDS,2)
%         for freq = 1:size(ERD_mean{1}{ch},2)
%             Dis_Eu_(freq,ch,:) = Dis_ERD_mean{ch}(freq,:);
%         end
%     end
%     addpath(genpath('G:\Dropbox\ERD\Codes\Topoplots'))
%     load HeadModel.mat
%     load G:\Dropbox\ERD\Codes\Topoplots\BCICIV_2a\electrodesBCICIV2a.mat
%     t_e    = 70;             % tamaño visual de los electrodos seleccionados.
%     M1.xy  = elec_pos(:,1:2);% posicion de los canales.
%     M1.lab = Channels;       % nombre de los canales.
%     type = 4;                % tipo de topoplot.
%     tex = 0;                 % 1 si quieres nombres de los canales.
%     sel = 1:22;              % Canales seleccionados.
%     database = 2;            % Base de datos seleccionada 2 es para BCICIV_2a
%     t_ = len*fs; %------------ time segment
%     ovlpt = round(0.9*t_);
%     ts = 1:t_-ovlpt:(7*fs)-t_;
%     set(0,'DefaultFigureWindowStyle','docked')
%     warning off
%     for fre = 1:5
%     for v = [1,5,10,14,18,22,25]
%         figure
%         rel = Dis_Eu_(fre,1:22,v)./max(max(max(Dis_Eu_)));
%         MyTopo_fun(rel,sel,M1.xy,M1.lab,[0 1],0,0,t_e,database,HeadModel,type,tex,colormap('parula'))
%         axis off square
%         colorbar
%         title(['Time: ' num2str(ts(v)/fs) ' - ' num2str(ts(v)/fs+len), ' Max: ' num2str(max(rel))])
%     end
%     end
% end


% close
%
% a = 1;
% fs = 250;
% t = (1:1750)/fs;
% % for ch = [8,12]
% %     figure
% %     SS = [2,8];
% %     for ss = SS
% %         for clas =1:2
% %             load(['F:\BCI\BCICIV_2a_0' num2str(ss) filesep  num2str(ss) '_' num2str(clas) '_'  'ERD_g.mat'],'erd')
% %             subplot(2,2,a)
% %             load erdscolormap
% % colormap(erdcolormap)
% %             ERD{clas}(:,:,ch) = erd.ERDS{ch}.erds';
% %             imagesc(t,1:17,squeeze(ERD{clas}(:,:,ch))*100,[-100,150])
% %             axis xy
% %             title(['Subject: ',num2str(ss),'Ch: ',num2str(ch)])
% %             a = a+1;
% %             if ss == 8
% %                 xlabel(['Class: ',num2str(clas)])
% %             end
% %         end
% %     end
% %     a = 1;suptitle('ERD')
% % end
%
% posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % posiciones de las graficas.
% set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% SS = 2;
% load('F:\BCI\BCICIV_2a_08\8_1_5ERDfolds.mat')
% for fold = 1:5
%     figure
%     for s = SS
%         for ch = 1:22
%             subplot(6,7,posi(ch))
%             load erdscolormap
%             colormap(erdcolormap)
%             imagesc(erd{fold}.ERDS{ch}.erds', [-1, 1.5]);
%             if ch==1; aa{fold} = erd{fold}.ERDS{ch}.erds'; end
%         end
%     end
% end

% %% Average/variance
% r2 = calcAveVar(s, h, [0, 0.05, 8], 'class', classes, 'heading', name, 'montage', [1 1 1 1], 'cue', 3);
% plotAveVar(r2);
%
% %% Combination of ERDS maps and average/variance
% r3 = calcCombiMap(s, h, [0, 0.05, 8], [5, 18, 40], 'method', method, 'class', classes, 'ref', [0.5, 1.5], 'f_bandwidths', [2, 4], 'f_steps', [1, 1], 'sig', 'boxcox', 'lambda', 1, 'alpha', 0.05, 'heading', name, 'montage', [1 1 1 1], 'cue', 3);
% plotCombiMap(r3);