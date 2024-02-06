SS = [7];
time = 1:1/512:6;
fs = 512;
mmax =  1.5;%max([ERDsfilt_{fold}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
mmin  = -1;%min ([ERDsfilt_{fol
% d}{1}{8}(3,:),ERDsfilt_{fold}{1}{12}(3,:),ERDsfilt_{fold}{2}{8}(3,:),ERDsfilt_{fold}{2}{12}(3,:)]);
for s = SS
    load(['I:\ERD_giga\MI\30_folds_training\ERD_folds_s',num2str(s),'.mat'])
    for freq = [3]
        a = 1;
        for ch = [13,48,50,21,31,58]
%             figure(a)
            set(gcf,'position',[560   528   560   420])
            x1 = 0.19;     y1 = 0.1827;     w = 0.775;      h = 0.7423;
            subplot('position',[x1,y1,w,h]);
            
            a = a+1;
            chan = ch;
            inBetween = [[-1.2, 1.6], [1.6, -1.2]];
            fill([[-0.1, -0.1],[2.5, 2.5]], inBetween, [244/255,244/255,244/255],'EdgeColor',[244/255,244/255,244/255]);
            hold on  
            ERD_{1}{ch} = mean(ERDsfilt_o{1}{1}{chan}([3,5,7,9],1*fs:6*fs),1);
            ERD_{2}{ch} = mean(ERDsfilt_o{1}{2}{chan}([3,5,7,9],1*fs:6*fs),1);
            %             subplot(2,3,1)
            %             plot(time,ERDsfilt_o{1}{1}{chan}(freq,1*fs:6*fs),'r','LineWidth',1.8)
            %             plot(time,ERDsfilt_o{1}{1}{chan}(freq,1*fs:6*fs),'Color',[255 180 128]/255,'LineWidth',1.8,'LineStyle','--')
            plot(time,ERD_{1}{ch},'Color',[0.8500 0.3250 0.0980],'LineWidth',1.8,'LineStyle','-')
            
            %             plot(time,ERDsfilt_o{1}{2}{chan}(freq,1*fs:6*fs),'Color',[0 153 0]/255,'LineWidth',1.8)
            %             plot(time,ERDsfilt_o{1}{2}{chan}(freq,1*fs:6*fs),'Color',[136 96 208]/255,'LineWidth',1.3)
            plot(time,ERD_{2}{ch},'Color',[0.9290 0.6940 0.1250],'LineWidth',1.3)
            %             title(['Freq ',num2str(freq),' chan:C3'])
            ylim([mmin mmax])
            xlim([1,6])
            set(gca,'TickLabelInterpreter','latex','Fontsize', 20)
            if s == 22
                xlabel('Time [$\it S$]' ,'Interpreter','latex')
            end
            if ch == 21 && s==14
                ylabel('ERD/ERS[%]','Interpreter','latex','Position',[0.15 0.25 -1])
            end
            if s== 22 && ch  ==21
                ylabel('ERD/ERS[%]','Interpreter','latex','Position',[0.15 0.25 -1])
            end
%             title(['channel ',num2str(ch)])
            %             saveas(gca,['Subject ',num2str(s),' Freq ',num2str(freq),' Ch ',num2str(ch)],'png')
%             saveas(gca,['F:\Dropbox\[1a] ERD Session\Figure\ERDs_plots_giga',filesep,'ERD_Sub_',num2str(s),'_ch_',num2str(ch)],'epsc')
%             saveas(gca,['F:\Dropbox\[1a] ERD Session\Figure\ERDs_plots_giga',filesep,'ERD_Sub_',num2str(s),'_ch_',num2str(ch)],'png')
%             close
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