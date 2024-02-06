SS = [14,22];
time = 1:1/512:6;
fs = 512;
mmax =  1.5;%max([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
mmin  = -1;%min ([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
for s = SS
    load(['L:\ERD_giga\MI\30_folds_training',filesep,'ERD_folds_s',num2str(s),'.mat'])
    for freq = [3,7,9]
        for ch = [13,48,50,21,31,58]
            
            figure
            chan = ch;
            %             subplot(2,3,1)
%             plot(time,ERDsfilt_o{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            plot(time,ERDsfilt_o{1}{1}{chan}(freq,1*fs:6*fs),'Color',[255 180 128]/255,'LineWidth',1.8,'LineStyle','--')
            hold on
%             plot(time,ERDsfilt_o{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            plot(time,ERDsfilt_o{1}{2}{chan}(freq,1*fs:6*fs),'Color',[136 96 208]/255,'LineWidth',1.3)
            %             title(['Freq ',num2str(freq),' chan:C3'])
            ylim([mmin mmax])
            xlim([1,6])
            set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            saveas(gca,['Subject ',num2str(s),' Freq ',num2str(freq),' Ch ',num2str(ch)],'png')
            saveas(gca,['F:\Dropbox\[1a] ERD Session\Figure\ERDs_plots_giga',filesep,'ERD_Sub_',num2str(s),'_ch_',num2str(ch)],'epsc')
%             close
            %         %     figure
            %         chan = 48;
            %         subplot(2,3,2)
            %         plot(time,ERDsfilt_{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %         hold on
            %         plot(time,ERDsfilt_{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %         title(['Freq ',num2str(freq),' chan:Cz'])
            %         xlim([1,6])
            %         set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            %
            %         %     figure
            %         chan = 50;
            %         subplot(2,3,3)
            %         plot(time,ERDsfilt_{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %         hold on
            %         plot(time,ERDsfilt_{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %         title(['Freq ',num2str(freq),' chan:C4'])
            %         xlim([1,6])
            %         set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            %
            %         %     figure
            %         chan = 21;
            %         subplot(2,3,4)
            %         plot(time,ERDsfilt_{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %         hold on
            %         plot(time,ERDsfilt_{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %         title(['Freq ',num2str(freq),' chan:P3'])
            %         xlim([1,6])
            %         set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            %
            %         %     figure
            %         chan = 31;
            %         subplot(2,3,5)
            %         plot(time,ERDsfilt_{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %         hold on
            %         plot(time,ERDsfilt_{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %         title(['Freq ',num2str(freq),' chan:Pz'])
            %         xlim([1,6])
            %         set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            %
            %         %     figure
            %         chan = 58;
            %         subplot(2,3,6)
            %         plot(time,ERDsfilt_{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %         hold on
            %         plot(time,ERDsfilt_{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %         title(['Freq ',num2str(freq),' chan:P4'])
            %         xlim([1,6])
            %         set(gca,'TickLabelInterpreter','latex','Fontsize', 15)
            % %         saveas(gca,['Subject ',num2str(s),' Freq ',num2str(freq)],'fig')
            %         saveas(gca,['Subject ',num2str(s),' Freq ',num2str(freq)],'epsc')
            %         saveas(gca,['F:\Dropbox\ERD\ERD_',filesep,'ERD_Sub_',num2str(s),'_ch_'])
            %             suptitle(['Subject ',num2str(s)])
        end
    end
end