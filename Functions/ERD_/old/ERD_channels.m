set(0,'DefaultFigureWindowStyle','docked')
SS = [8,1,7];
for s = SS
    % lf
    % load(['D:\BCI\ERDSjunio2019\ERD_folds10sin_sub' num2str(s) '.mat'])
    % fy
    %load(['D:\BCI\ERDS sin significancia\ERD_folds10sin_sub' num2str(s) '.mat'])
    load(['E:\BCI\ERDSjunio2019\ERD_folds10sin_sub' num2str(s) '.mat'])
    fold = 1;
    time = 0:1/250:7;
    mmax =  0.65;%max([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
    mmin  = -0.72;%min ([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
    
    for chan = [8,12]
        figure;
        %set(gcf,'position',[667   528   600   300])
        %x1 = 0.15;     y1 = 0.1;     w = 0.8;      h = 0.85;
        %subplot('position',[x1,y1,w,h]);
        
        plot(time,ERDsfilt_{fold}{1}{chan}(3,:),'r','LineWidth',1.8)
        hold on
        plot(time,ERDsfilt_{fold}{2}{chan}(3,:),'Color',[0 153 0]/255,'LineWidth',1.8)
        set(gca,'TickLabelInterpreter','latex')
        set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
        ylim([mmin-0.1 mmax])%mmax+0.1
        xlim([1,6])
        saveas(gca,['ERD_Sub_' num2str(s) '_ch_' num2str(chan) ],'epsc')
%         saveas(gca,['Sub_' num2str(s) '_class_' num2str(chan) ],'png')
%         close
    end
end