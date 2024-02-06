set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [1 2];                                                               % posicion de grafica.
% Rango de freq.
frang = [12 25];                                                             % rango de frecuencias.
% Rango de tiempo
trang = [16 116];                                                           % rango de tiempo.
% rango de graficar limites
rt = [1 5]; % tiempo                                                        % rango de tiempos a graficar.
rf = [6 12];% freq                                                          % rango de frecuencias a graficar.
% mejor sujeto
s = [3];                                                                    % mejor sujeto de la base de datos.
% canal Cz
ChanCz = 10;
% graficador en 3D todos los canales
for sub = 1:9 % numero de la figura
    % grafica 1
    if sub == 3
        sub = 3;
    else
        figure(sub+20)
        subplot(1,2,posi(1))
        surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{s(1)}(ChanCz,frang(1):frang(2),trang(1):trang(2))))
        title(['ch: ' labels{ChanCz} ' Sujeto: ' num2str(s(1))])
        % title(['ch: ' num2str(ChanCz) ' Sujeto: ' num2str(s(1))])
        xlim([rt(1) rt(2)])
        ylim([rf(1) rf(2)])
        shading interp
        view(-18,77)
        colormap jet
        xlabel(['time_{seconds}'])
        ylabel('Freq')
        zlabel('Amp')
        hold on
        plot3(t.*0+2,t.*5,t.*0,'r','LineWidth',2)
        plot3(t.*0+4.5,t.*5,t.*0,'r','LineWidth',2)
        hold off
        drawnow
        
        % grafica 2
        subplot(1,2,posi(2))
        surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(ChanCz,frang(1):frang(2),trang(1):trang(2))))
        title(['ch: ' labels{ChanCz} ' Sujeto: ' num2str(sub)])
        % title(['ch: ' num2str(ChanCz) ' Sujeto: ' num2str(s(2))])
        xlim([rt(1) rt(2)])
        ylim([rf(1) rf(2)])
        shading interp
        view(-18,77)
        colormap jet
        xlabel(['time_{seconds}'])
        ylabel('Freq')
        zlabel('Amp')
        hold on
        plot3(t.*0+2,t.*5,t.*0,'r','LineWidth',2)
        plot3(t.*0+4.5,t.*5,t.*0,'r','LineWidth',2)
        hold off
        drawnow
    end
end
