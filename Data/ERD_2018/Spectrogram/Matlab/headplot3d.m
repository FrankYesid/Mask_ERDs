load('ERD.mat')
%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 3D %%%%%%%%%%% %%% con imagesc or surf
load('posiciones.mat')
% M1.lab = {'FP1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','Poz','Pz','CPz','FPz','FP2','AF8','AF4','Afz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'}';
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
% posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
% % ubicación de cada una de las graficas. para BCICIB_2a.mat
% Rango de freq.
% frang = [5 40];                                                             % rango de frecuencias a graficar.
% Rango de tiempo
% trang = [5 595];                                                           % rango de tiempo a graficar.
% rango de graficar limites
rt = [1 6]; % tiempo                                                        % rango en el plot del tiempo.
rf = [4 40];% freq                                                          % rango en el plot de las frecuencias.
pos = [5 15 16 27 26 25 24 35 36 37 38 49 48 47 46 57 58 59 60 71 70 69 68 67 81 82 93 105 94 83 72 61 6 7 19 18 17 28 29 30 31 32 43 42 41 40 39 50 51 52 53 54 65 64 63 62 73 74 75 76 77 85 84 95];
% graficador en 3D todos los canales
for sub = [3,8,14,16,17,38,40,41,43,50]
    i = 1;
    fig = figure(sub);
    for ch = 1:N_channels
        subplot(10,11,pos(i))
        %         surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(ch,frang(1):frang(2),trang(1):trang(2))))
        surf(t,f,squeeze(ERD{sub}(ch,:,:)))
        %         imagesc(t,f,squeeze(ERD{sub}(ch,:,:)))
        %         title([labels{ch}])
        title(M1.lab(ch),'Interpreter','latex')
        xlim([rt(1) rt(2)])
        ylim([rf(1) rf(2)])
        yticks([])
        
        shading interp
        %         view(-18,38)
        view(0,90)
        colormap jet
        grid on
        hold on
        line([2.5,2.5],[1,max(f)],'LineWidth',2,'color','r')
        line([4.5,4.5],[1,max(f)],'LineWidth',2,'color','r')
        hold off
        if pos(ch) == 105
            xlabel('time (s)','Interpreter','latex')
            ylabel('Freq','Interpreter','latex')
            zlabel('Amp','Interpreter','latex')
            yticks([4,40])
        end
        drawnow
        i = i+1;
    end 
    suptitle(['Subjet '  num2str(sub)])
    saveas(fig,['Subject ' num2str(sub) '.png'])
    saveas(fig,['Subject ' num2str(sub) '.fig'])
% annotation('textbox', [0 0.9 1 0.1],'String', 'hello, title','EdgeColor', 'none','HorizontalAlignment','center','Interpreter','latex')
    %     grafica de label(x,y,z)
    %     subplot(6,7,42)
    %     surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(10,frang(1):frang(2),trang(1):trang(2)).*0))
    %     xlabel('time (s)','Interpreter','latex')
    %     ylabel('Freq','Interpreter','latex')
    %     zlabel('Amp','Interpreter','latex')
    %     shading interp
    %     view(-18,38)
    %     xlim([rt(1) rt(2)])
    %     ylim([rf(1) rf(2)])
    %     grid on
    %     hold on
    %     plot3(t.*0+2,t.*5,t.*0,'-r','LineWidth',1.5)
    %     plot3(t.*0+4.5,t.*5,t.*0,'-r','LineWidth',1.5)
    %     hold off
    %     colormap jet
end