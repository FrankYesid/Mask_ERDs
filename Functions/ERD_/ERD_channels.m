% set(0,'DefaultFigureWindowStyle','docked')
SS = [8];
con = 1;
freq= 3;
freq1= [9];
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
for s = SS
    % lf
    % load(['D:\BCI\ERDSjunio2019\ERD_folds10sin_sub' num2str(s) '.mat'])
    % fy
    %load(['D:\BCI\ERDS sin significancia\ERD_folds10sin_sub' num2str(s) '.mat'])
    %     load(['E:\BCI\ERDSjunio2019\ERD_folds10sin_sub' num2str(s) '.mat'])
    if con == 1
        load(['I:\ERD\all_trials\ERD_s',num2str(s),'.mat'],'ERDsfilt_')
    else
        load(['I:\ERD\all_trials\ERD_s',num2str(s),'_sin.mat'],'ERDsfilt_')
    end
    %     fold = 1;
    time = 0:1/250:7;
    mmax =  0.65;%max([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
    mmin = -0.72;%min ([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
    for chan = [8:12]
        figure;
        set(gcf,'position',[667   128   600   300])
        x1 = 0.15;  y1 = 0.1;     w = 0.8;      h = 0.85;
        subplot('position',[x1,y1,w,h]);
        %         subplot(1,5,chan-7)
        %         plot(time,ERDsfilt_{1}{chan}(:,freq),'r','LineWidth',1.8)
        hold on
        plot(time,ERDsfilt_{1}{chan}(:,freq),'Color',[0.3,0.75,0.93],'LineWidth',1.8)
        plot(time,squeeze(mean(ERDsfilt_{1}{chan}(:,freq1),2)),'Color',[0.3,0.75,0.93],'LineWidth',1.8,'LineStyle','--')
        %         erd_p = ERDsfilt_{1}{chan}(:,freq)';
        %         patch([time nan],[erd_p nan],[ones(1,1751) nan],[erd_p nan],...
        %             'FaceColor','none','EdgeColor','interp','LineWidth',1.8)
        set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
        ylim([-1 1.5])
        %         caxis([-1 1.5])
        %         colormap(erdcolormap)
        %         subplot(2,5,chan-2)
        %         erd_p2 = ERDsfilt_{2}{chan}(:,freq)';
        %         patch([time nan],[erd_p2 nan],[ones(1,1751) nan],[erd_p2 nan],...
        %             'FaceColor','none','EdgeColor','interp','LineWidth',1.8)
        %         plot(time,ERDsfilt_{2}{chan}(:,freq),'Color',[0 153 0]/255,'LineWidth',1.8)
        plot(time,ERDsfilt_{2}{chan}(:,freq),'Color',[0,0.45,0.74],'LineWidth',1.8)
        plot(time,squeeze(mean(ERDsfilt_{2}{chan}(:,freq1),2)),'Color',[0,0.45,0.74],'LineWidth',1.8,'LineStyle','--')
        set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
        ylim([-1 1.5])
        line([0,7],[0,0],'LineStyle','--','Color',[0.8,0.8,0.8],'LineWidth',1.8)
        line([2,2],[-1,1.5],'LineStyle','--','Color',[0.8471,0.2157,0.1725],'LineWidth',1.8)
        %         caxis([-1 1.5])
        %         colormap(erdcolormap)
        
        %         if chan == 8 || chan == 9 || chan == 10
        %             text(3,max(erd_p),'Clase 1','Interpreter','latex','Fontsize', 15)
        %             text(3,min(erd_p2),'Clase 2','Interpreter','latex','Fontsize', 15)
        %         elseif chan == 11 || chan == 12
        %             text(3,min(erd_p),'Clase 1','Interpreter','latex','Fontsize', 15)
        %             text(3,max(erd_p2),'Clase 2','Interpreter','latex','Fontsize', 15)
        %         end
        
        % ylim([mminv0.1 mmax])%mmax+0.1
        %         % xlim([1,6])
%         if con == 1
%             saveas(gca,['ERD_Sub_' num2str(s) '_ch_' num2str(chan) '_Con_lab'],'epsc')
%         else
%             saveas(gca,['ERD_Sub_' num2str(s) '_ch_' num2str(chan) '_sin_lab'],'epsc')
%         end
        % saveas(gca,['Sub_' num2str(s) '_class_o' num2str(chan) ],'png')
        %matlab2tikz(['C:\Users\frany\Desktop\Documento_tesis\Figures\Figs_preliminar\con_lap_ordenado\curve_ERD_S_',num2str(s),'_ch',num2str(chan),'.tikz'])
        %close
    end
end

% organizar una versión mas amplia de todos los canales, colocando en otra
% vista para ver la gráfica en version 3D.
