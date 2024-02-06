%% Programa ERD Spectrogram & Wavelet
%-------------------------------------------------------------------------
% Inputs:
% X{subject}{trial}(samples x chan): Data to analysis
% y{subject}{trial}
% fs: Sample frequency
% labels{}
%-------------------------------------------------------------------------
% Outputs:
% Xf: Filtered data
% Xf{trial}(samples x chan)
% segun el uso de spectrogram o wavelet packet (con maximo o con entropia).
% ERD{subject}(chan x freq x time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 Signal Processing and Recognition Group
% F.Y. Zapata C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seccion A
%% %%%%%%%%%%%%%%%%%%%%%%%% ERDs - SPECTROGRAM %%%%%%%%%%%%%%%%%%%%%%%% %%
%% Limpiar
clear all; clc                                  % limpia la ventada de comandos y también la memoria del espacio de trabajo.
%% Load Database
load('BCICIV_2a.mat');                          % carga base de datos EEG.
load('labels.mat');                             % carga base de datos de los canales mas relevantes.
load('layout');                                 % carga base de datos.
%% Filter Band
clc;
a = [13 30];                                    % rango del filtro.
A = cell(9,1);                                  % celdas de la señal filtrada.
% filtrado de cad sujeto.
for k = 1:9
    fprintf(['Sujeto: ' num2str(k) ' de 9' '\n'])
    Xfreq = fcnfiltband(X{k},fs,a,5);           % funcion Filtro.
    A{k} = Xfreq;
end
X = A;                                          % almacena la señal filtrada.
%% Prueba de ERD con diferentes bandas de R e {2-40} Hz con traslape de 1 Hz.
%% Filter Band
cfil = cell(38,1);
%%
% for fil = 11:15
a = [4 40];                          % rango del filtro.
%     A = cell(9,1);                              % celdas de la señal filtrada.
% filtrado de cad sujeto.
for k = 1:9
    %         fprintf(['Filtro: ' num2str(fil+1) ' a ' num2str(fil+3) ' Hz' ' Sujeto: ' num2str(k) ' de 9' '\n'])
    Xfreq{k} = fcnfiltband(X{k},fs,a,5);       % funcion Filtro.
    %         A{k} = Xfreq;
end
%     cfil{fil} = A;
%     save(['matriz' num2str(fil) '.mat'],'A')
%     clear A
% end
clear X
X = Xfreq;
clear Xfreq
%% %%%%%%%%%%%%%%% Permutation entropy (fast algorithm) %%%%%%%%%%%%%%% %%
% windowSize = 0.5*fs;                          % ventana.
tau = 1;                                        % - time lag scalar OR time lag vector (length = n-1)
order = 3;                                      % - permutation order
% acc = ;
METHOD = 'noise';                               % method - method how to treat equal values
PerEnt = cell(9,1);                             % celdas almacenadoras del tiempo segun la entropia.
HEntro = cell(9,1);                             % celdas almacenadoras de entropia de cada canal.
Clase = [1];                                    % clase o clases seleccionadas.
tt1 = 250;                                      % tiempo de referencia.
Canales = [1:22]; N_canal = length(Canales);    % cantidad de canales.
% ciclo para la entropia
for sub = 1:9
    data = X{sub}(y{sub} == Clase);                                         % informacion sobre cada clase (trials segun la clase).
    N_trial = numel(data);                                                  % cantidad de trials.
    PerEnt{sub} = zeros(N_trial,N_canal,1);                                 % tamaño de la matriz de permutacion entropy.
    HEntro{sub} = zeros(N_trial,N_canal,1);                                 % tamaño de la matriz de permutacion entropy en cada channel.
    for tri = 1:N_trial
        fprintf(['sujeto: ' num2str(sub) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        for cnl = 1:N_canal
            % signal
            signal = data{tri}(tt1:end,Canales(cnl));                       % señal de 1 a 5 segundos.
            % Calcular Entropy
            [HE,mat1,mat2] = PermuEntropy(signal,order,tau,METHOD);
            % mat = [mat1 mat2];
            dmax = max(mat2);
            dmax = find(mat2 == dmax);
            dmax = max(dmax);
            Hentch = mat1(dmax);
            % Almacenar entropias
            HEntro{sub}(tri,cnl,:) = HE;
            PerEnt{sub}(tri,cnl,:) = Hentch+tt1;                            % Trials x canal x instante
        end
    end
end
% promedio de los tiempos de referencia.
fprintf(['En{P_entropy(t,f)}\n'])
for sub = 1:9
    fprintf(['sujeto: ' num2str(sub) ' de ' '9' '\n'])
    % En{Pnc(t,f)}
    Pent_c{sub,1} = (squeeze(mean(PerEnt{sub},1)))./250;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% ERDs spectrogram %%%%%%%%%%%%%%%%%%%%%%%%% %% % parametros de la STFT
clc;
% Nw = 1;                                       % tama?o de la ventana en seg
Nw = fs;                                        % tamaño de la ventana.
window = hamming(Nw);                           % Vector de la ventana.
noverlap = floor(0.98*Nw);                      % interpolación de la ventana.
X_suj = cell(9,1);                              % celdas almacenadoras de spectrogram.
Clase = [2];                               % clase o clases seleccionadas.
for su = 1:9                                % calculo del espectrograma relacionada a EEG.
    Sujeto = X{su};                                                         % sujeto seleccionado.
    etiqueta = y{su};                                                       % etiqueta seleccionada.
    % selecciona los trail de las etiquetas solicitadas.
    % Izq_Der = Sujeto(ismember(etiqueta,[1 2]));
    Izq_Der = Sujeto(ismember(etiqueta,Clase));                             % seleccion de los trails segun la clase o clases seleccionadas.
    % escojer el que tenga la menor cantidad de trials
    N_trial = length(Izq_Der);                                              % cantidad de trials.
    Canales = [1:22]; N_canal = length(Canales);                            % cantidad de canales.
    % Inicializar matrices
    %     temp = spectrogram(Izq_Der{1}(16:end,1), window, noverlap,Nw);               % tamaño temporal de los datos en cada sujeto.
    %     temp = size(temp,2);
    %     X_suj(su) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)};              % Trials x canal x frecuencia x tiempo.
    for tri = 1:N_trial
        fprintf(['sujeto: ' num2str(su) ' de 9' ' Class: ' num2str(Clase) ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        xmean = mean(mean(Izq_Der{tri}));
        for cnl = 1:N_canal
            % signal
            signal = (Izq_Der{tri}(:,Canales(cnl)));
            %             signal = signal-xmean;
            % Calcular STFT
            [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
            
            X_Class = abs(X_Class);                                         % absoluto para tener la señal (freq,tiempo)
            
            % Almacenar espectrogramas
            X_suj{su}(tri,cnl,:,:) = X_Class;                               % Trials x canal x frecuencia x tiempo
        end
    end
end
% X_1 = X_suj;
% clearvars -EXCEPT X_1 t f
%%
for su = [3,8]
    for tr = 1:numel(X_suj{su})
        for fil = 1:size(X_suj{su}{tr},2)
            for ch = 1:size(X_suj{su}{tr},1)
                Dat{su}{fil}{tr}(ch,:) = X_suj{su}{tclearr}(ch,fil,:);
            end
        end
    end
end
% tiempo de referencia.
% Sujetos = X_suj;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% ERD - Normal %%%%%%%%%%%%%%%%%%%%%%%%%%% %% % calculo del ERD
t1 = 0;                                         % tiempo inicial de referencia.
t2 = 2;                                         % tiempo final de referencia.
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - t2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
%%
ul = 300;
up = 500;

% r_nc = cell(9,1);                               % celda para esperanza respecto al tiempo de referencia.
% r_c = cell(9,1);                                % celda para esperanza respecto a los trials (intentos).
% m_c = cell(9,1);                                % celda para esperanza respecto a los trials (intentos).
%% mean ERD
su = 1:9;
class = Clase;
fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
for y =su
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % Et{Pnc(t,f):te[Ta,Tb]}
    r_nc{y,1} = squeeze(mean(X_suj{y}(:,:,:,ul:up),4));
end
fprintf(['En{r_nc(f)}\n'])
for y =su
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{r_nc(f)}
    r_c{y,1} = squeeze(mean(r_nc{y,1},1));
end
fprintf(['En{Pnc(t,f)}\n'])
for y =su
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{Pnc(t,f)}
    m_c{y,1} = squeeze(mean(X_suj{y},1));
end
fprintf(['ERD = m_c/r_c\n'])
for y = su
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}{class}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1; % channel,freq,time
end

% % c = 22;
% % figure(23)
% % plot(t,squeeze(ERD{sub}(c,8,:)))
% % save('ERD_t13.mat','ERD')
%%
clearvars -EXCEPT ERD t f
%% mean ERD - trial

for su = [3,8]
    for fil = 1:numel(Dat{su})
        for tr = 1:numel(Dat{su}{fil})
            r_nc{su}{fil}{tr}(ch) = squeeze(Dat{su}{fil}{ch,t1:t2});
            ERD{su}{fil}{tr}(ch,:) = bsxfun(@times,Dat{su}{fil}{tr}(ch,:),1./r_nc{su}{fil}{tr}(ch)) - 1;
        end
    end
end

%% ERD hilbert
Clase = [1];                                      % clase o clases seleccionadas.
for su = [3,2]                                   % calculo del espectrograma relacionada a EEG.
    Sujeto = Xfreq{su};                                                         % sujeto seleccionado.
    etiqueta = y{su};                                                       % etiqueta seleccionada.
    % selecciona los trail de las etiquetas solicitadas.
    % Izq_Der = Sujeto(ismember(etiqueta,[1 2]));
    Izq_Der = Sujeto(ismember(etiqueta,Clase));                             % seleccion de los trails segun la clase o clases seleccionadas.
    % escojer el que tenga la menor cantidad de trials
    N_trial = length(Izq_Der);                                              % cantidad de trials.
    Canales = [1:22]; N_canal = length(Canales);                            % cantidad de canales.
    % Inicializar matrices
    %     temp = spectrogram(Izq_Der{1}(16:end,1), window, noverlap,Nw);               % tamaño temporal de los datos en cada sujeto.
    %     temp = size(temp,2);
    %     X_suj(su) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)};              % Trials x canal x frecuencia x tiempo.
    for tri = 1:N_trial
        fprintf(['sujeto: ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        for cnl = 1:N_canal
            % signal
            signal = hilbert(Izq_Der{tri}(:,Canales(cnl)));
            
            % Calcular STFT
            %             [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
            
            X_Class = real(abs(signal));                                         % absoluto para tener la señal (freq,tiempo)
            
            % Almacenar espectrogramas
            X_suj{su}(tri,cnl,:) = X_Class;                               % Trials x canal x frecuencia x tiempo
        end
    end
end
%% %%%%%%%%%%%%%%%% prueba - ERD - señal de cada trial %%%%%%%%%%%%%%%% %% ---------- para base de datos (DATABASE)
t1 = 0;                                         % tiempo inicial de referencia.
t2 = 2;                                         % tiempo final de referencia.
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - 2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
r_c = cell(22,1);                               % celda para esperanza respecto a los trials (intentos).
m_c = cell(22,1);                               % celda para esperanza respecto a los trials (intentos).
ff = cell(9,1);
ERD = cell(9,1);
for sub = 1:numel(X_suj)
    tr = size(X_suj{sub},1);
    ERD{sub} = cell(tr,1);
    for trial = 1:tr
        fprintf(['sujeto: ' num2str(sub) ' de ' '9' ' ...trial: ' num2str(trial) ' de ' num2str(tr) '\n'])
        %         fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
        for y = 1:22
            %             fprintf(['Channel: ' num2str(y) ' de ' '22' '\n'])
            % Et{Pnc(t,f):te[Ta,Tb]}
            r_c{y,1} = squeeze(mean(X_suj{sub}(trial,y,:,ul:up),4));
        end
        
        %         fprintf(['En{Pnc(t,f)}\n'])
        for y = 1:22
            %             fprintf(['Channel: ' num2str(y) ' de ' '22' '\n'])
            % En{Pnc(t,f)}
            m_c{y,1} = squeeze(X_suj{sub}(trial,y,:,:));
        end
        
        %         fprintf(['ERD = m_c/r_c\n'])
        for y = 1:22
            %             fprintf(['Channel: ' num2str(y) ' de ' '22' '\n'])
            % ERD = m_c/r_c
            ERD_temp{y,1} = bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1;        % channel,freq,time
        end
        ERD_ch = cell(22,1);
        for y = 1:22
            ET = sum(ERD_temp{y,1}(9:21,16:116).^2,2);                      % selecciona freq(8-13) - tiempos(1-5)
            % señecciona la mejor potencia de la señal.
            dd = max(ET);
            dat = find(ET == dd) + 6;
            ff{sub}(trial,y,:) = f(dat);
            ERD_ch{y} = ERD_temp{y,1}(dat,:);
        end
        ERD_ch = cell2mat(ERD_ch);
        ERD{sub}{trial} = ERD_ch;
    end
end
%% Variables
ERDs = cell(2,1);
ffe = cell(2,1);
%% clase 1
ERDs{1} = ERD;
ffe{1} = ff;
%% clase 2
ERDs{2} = ERD;
ffe{2} = ff;
%%
%% guardas database
save('ERD','ERDs','t','ffe')                    % guarda la informacion de los datos de ERDs calculado.
%% %%%%%%%%%%%%%%%  ERD calculado referente a la entropy %%%%%%%%%%%%%%% %%
r_nc = cell(9,1);                               % celda para esperanza respecto al tiempo de referencia.
r_c = cell(9,1);                                % celda para esperanza respecto a los trials (intentos).
m_c = cell(9,1);                                % celda para esperanza respecto a los trials (intentos).
fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])           % calculo del ERD
for y = 1:9
    for ch = 1:22
        fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
        % tiempo de referencia.
        t1 = 0;
        t2 = Pent_c{y}(1,ch);
        temp1 = abs(t - t1);
        min1 = min(temp1);
        temp2 = abs(t - 2);
        min2 = min(temp2);
        ul = find(temp1 == min1);
        up = find(temp2 == min2);
        % calculo de ERD
        r_nc{y,1} = squeeze(mean(X_suj{y}(:,ch,:,ul:up),4));
    end
end
fprintf(['En{r_nc(f)}\n'])
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{r_nc(f)}
    r_c{y,1} = squeeze(mean(r_nc{y,1},1));
end
fprintf(['En{Pnc(t,f)}\n'])
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{Pnc(t,f)}
    m_c{y,1} = squeeze(mean(X_suj{y},1));
end
fprintf(['ERD = m_c/r_c\n'])                    % calculo del ERD.
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1; % channel,freq,time
end

%% Marginal de ERD
merd = cellfun(@(x) squeeze(mean(x(:,:,:),2)),ERD,'UniformOutput',false);

%% Plot de Marginal - ERD
ns = 8;
cha1 = 8;
cha2 = 12;

figure('name','ERDs')
subplot(2,1,1)
plot(t,merd{ns}(cha1,:)), title('Channel C3'),
hold on, plot([2 2],[min(merd{ns}(cha1,:)) max(merd{ns}(cha1,:))]), plot([3.25 3.25],[min(merd{ns}(cha1,:)) max(merd{ns}(cha1,:))]),plot([6 6],[min(merd{ns}(cha1,:)) max(merd{ns}(cha1,:))])
hold off
subplot(2,1,2)
plot(t,merd{ns}(cha2,:)), title('Channel C4'),
hold on, plot([2 2],[min(merd{ns}(cha2,:)) max(merd{ns}(cha2,:))]), plot([3.25 3.25],[min(merd{ns}(cha2,:)) max(merd{ns}(cha2,:))]),plot([6 6],[min(merd{ns}(cha2,:)) max(merd{ns}(cha2,:))])
hold off

%% ERSP
ERSP = cellfun(@(X) squeeze(mean(X(:,:,:,:),1)),X_suj,'UniformOutput',false);

mersp = cellfun(@(X) squeeze(mean(X(:,:,:),2)),ERSP,'UniformOutput',false); % marginal ERSD
ns = 8;
cha1 = 8;
cha2 = 12;

figure('name','ERSP')
subplot(2,1,1)
plot(t,mersp{ns}(cha1,:)), title('Channel C3'),
hold on, plot([2 2],[min(mersp{ns}(cha1,:)) max(mersp{ns}(cha1,:))]), plot([3.25 3.25],[min(mersp{ns}(cha1,:)) max(mersp{ns}(cha1,:))]),plot([6 6],[min(mersp{ns}(cha1,:)) max(mersp{ns}(cha1,:))])
hold off
subplot(2,1,2)
plot(t,mersp{ns}(cha2,:)), title('Channel C4'),
hold on, plot([2 2],[min(mersp{ns}(cha2,:)) max(mersp{ns}(cha2,:))]), plot([3.25 3.25],[min(mersp{ns}(cha2,:)) max(mersp{ns}(cha2,:))]),plot([6 6],[min(mersp{ns}(cha2,:)) max(mersp{ns}(cha2,:))])
hold off

%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram (Firmas) - 2D %%%%%%%%%%% %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % posiciones de las graficas.
% rango de amplitud
Lrang = [-1.5 1.5];                                                         % rango de visualizacion en emplitud.
% graficador en 2D todos los canales
% freq = 14;
for sub = 1:9
    i = 1;
    figure
    for ch = 1:22
        subplot(6,7,posi(i))
        imagesc(t',f,squeeze(ERD{ch}),[Lrang(1) Lrang(2)])
        area(t(16:end),squeeze(ERD{sub}(ch,freq,16:end)))
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        xlim([1 5])
        %         ylim([1 40]);
        %         legend(num2str(freq))
        colormap jet
        drawnow
        % if s(sub,1) == 1
        %     title(['Ch: ' num2str(1)])
        % end
        i = i+1;
    end
end
%%
label_save = ['Canal_' num2str(canal)];
saveas(eval('fig1'),[label_save '.eps'])
saveas(eval('fig2'),[label_save '.eps'],'psc2');

%%  Analisis de relevancia -------------NO---------------------------------
% % load('X_sujetos_Izq_Der.mat')
% % Sujetos = f_pca(X_sujetos_class_1);
% % save('Sujetos.mat','Sujetos');
%% Graficar ---------------------------------------------------------------
% % clc
% % fs = 250;
% % Nw = 1; % tama?o de la ventana en seg
% % Nw = 1*fs;
% % window = hamming(Nw);
% % noverlap = floor(0.9*fs);
% % load('ERD.mat')
% % paso = Nw-noverlap;
% % load('Sujetos.mat')
% % [N_freq N_t] = size(squeeze(Sujetos{1}(1,1,:,:)));
% % f = [0:floor((N_freq-1))]*fs/(2*N_freq);
% % time = [0:N_t-1]*paso/fs+(Nw/2)/fs;
%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 2D %%%%%%%%%%% %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.
% rango de amplitud
% Lrang = [-2.5 2.5];                                                             % rango de colocar la señal.
% graficador en 2D todos los canales
% for sub = 1:9
% clas =1;
for sub = 3
    i = 1;
    for class = 1
        figu = figure;
        for ch = 1:22
            subplot(6,7,posi(ch))
            %         imagesc(t,f,squeeze(ERD{sub}(ch,1:21,95:595)),[Lrang(1) Lrang(2)])
            %             surf(t,f,(squeeze(ERD{sub}{class}(ch,:,:)))) % /max([max(max(max((ERD{sub}{1}(:,:,:))))), max(max(max((ERD{sub}{2}(:,:,:)))))])))s  %,[Lrang(1) Lrang(2)])
            %             surf(t,f,squeeze(X_72(ch,:,:)/max([max( max( max(X_71))), max( max( max(X_72)))])))
            %             errorbar(squeeze(mean(X_1{1}(:,ch,:,:),1)),squeeze(mean(X_2{1}(:,ch,:,:)),1))
            surf(t,f,squeeze(std(X_1{sub}(:,ch,:,:))))
            %             set(gca,'ColorScale','log')
            title(['Ch: ' num2str(ch)])
            shading interp
            %             set(gca,'ZScale','log')
            %             zlim([0 1])
            %             set(gca,'YScale','log')
            %             caxis([0 1])
            view([0 90])
            colorbar
            if ch == 22
                xlabel('time (seg)')
            end
            %         ylim([1 40])
            xlim([1 t(end)])
            ylim([1 30]);
            colormap jet
            drawnow
            % if s(sub,1) == 1
            %     title(['Ch: ' num2str(1)])
            % end
            i = i+1;
        end
        %         suptitle(['Sujeto ' num2str(sub) ' clase ' num2str(clas)])
        suptitle(['Sujeto ' num2str(sub) ' - clase ' num2str(class)])
        %         saveas(figu,['C:\Users\frany\Desktop\Frank\S' num2str(sub) 'c' num2str(class) 'std_.fig']);
        %         close all
    end
end

%% Headplot for channel ciclo (std) - potencial %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.

for sub = 3
    i = 1;
    for class = 1
        figu = figure;
        for ch = 1:22
            subplot(6,7,posi(ch))
            for freq = 5:41
                hold on
                errorbar(mean(squeeze(X_1{sub}(:,ch,freq,:)),1),std(squeeze(X_1{sub}(:,ch,freq,:)),1))
            end
            %             set(gca,'ColorScale','log')
            title(['Ch: ' num2str(ch)])
            shading interp
            %             set(gca,'ZScale','log')
            %             zlim([0 1])
            %             set(gca,'YScale','log')
            %             caxis([0 1])
%             view([0 90])
%             colorbar
            if ch == 22
                xlabel('time (seg)')
            end
            %         ylim([1 40])
%             xlim([1 t(end)])
%             ylim([1 30]);
            colormap jet
            drawnow
            % if s(sub,1) == 1
            %     title(['Ch: ' num2str(1)])
            % end
            i = i+1;
            hold off
        end
        %         suptitle(['Sujeto ' num2str(sub) ' clase ' num2str(clas)])
        suptitle(['Sujeto ' num2str(sub) ' - clase ' num2str(class)])
        %         saveas(figu,['C:\Users\frany\Desktop\Frank\S' num2str(sub) 'c' num2str(class) 'std_.fig']);
        %         close all
    end
end

%% %%%%%%%%%%%% Headplot for channel ciclo spectrogram - 3D %%%%%%%%%%% %%%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % ubicación de cada una de las graficas.
% Rango de freq.
frang = [5 16];                                                             % rango de frecuencias a graficar.
% Rango de tiempo
trang = [5 595];                                                           % rango de tiempo a graficar.
% rango de graficar limites
rt = [1 5]; % tiempo                                                        % rango en el plot del tiempo.
rf = [5 30];% freq                                                          % rango en el plot de las frecuencias.
% graficador en 3D todos los canales
for sub =3
    i = 1;
    for ch = 1:22
        figure
        subplot(6,7,posi(i))
        surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(ch,frang(1):frang(2),trang(1):trang(2))))
        %         title([labels{ch}])
        title(['ch: ' num2str(ch)])
        xlim([rt(1) rt(2)])
        ylim([rf(1) rf(2)])
        shading interp
        view(-18,38)
        colormap jet
        grid on
        hold on
        plot3(t.*0+2,t.*5,t.*0,'-r','LineWidth',1.5)
        plot3(t.*0+4.5,t.*5,t.*0,'-r','LineWidth',1.5)
        hold off
        drawnow
        i = i+1;
    end
    %     grafica de label(x,y,z)
    subplot(6,7,42)
    surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(10,frang(1):frang(2),trang(1):trang(2)).*0))
    xlabel(['time_{seconds}'])
    ylabel('Freq')
    zlabel('Amp')
    shading interp
    view(-18,38)
    xlim([rt(1) rt(2)])
    ylim([rf(1) rf(2)])
    grid on
    hold on
    plot3(t.*0+2,t.*5,t.*0,'-r','LineWidth',1.5)
    plot3(t.*0+4.5,t.*5,t.*0,'-r','LineWidth',1.5)
    hold off
    colormap jet
end
%% %%%%%%%%%%% Plot individual 3D - 2 sujetos(3,S) - canal CZ %%%%%%%%%% %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [1 2];                                                               % posicion de grafica.
% Rango de freq.
frang = [7 25];                                                             % rango de frecuencias.
% Rango de tiempo
trang = [16 116];                                                           % rango de tiempo.
% rango de graficar limites
rt = [1 5]; % tiempo                                                        % rango de tiempos a graficar.
rf = [1 30];% freq                                                          % rango de frecuencias a graficar.
% mejor sujeto
s = [3];                                                                    % mejor sujeto de la base de datos.
% canal Cz
ChanCz = 10;
load labels
% graficador en 3D todos los canales
for sub = 1:9 % numero de la figura
    % grafica 1
    if sub == 3
        sub = 3;
    else
        figure(sub+40)
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

%% %%%%%%%%%%% Plot individual 3D - todos sujetos - canal CZ %%%%%%%%%% %%
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% Rango de freq.
frang = [5 25];                                                             % rango de frecuencias.
% Rango de tiempo
trang = [16 116];                                                           % rango de tiempo.
% rango de graficar limites
rt = [1 5]; % tiempo                                                        % rango de tiempos a graficar.
rf = [5 30];% freq                                                          % mejor sujeto de la base de datos.
% canal Cz
ChanCz = 10;
% carga los nombres de cada punto de la cabeza.
load labels
% vista de la grafica.
Az = -14;
El = 77;
% graficador en 3D todos los canales
for sub = 1:9 % numero de la figura
    % grafica
    figure(sub)
    surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(ChanCz,frang(1):frang(2),trang(1):trang(2))))
    title(['ch: ' labels{ChanCz} ' Sujeto: ' num2str(sub)])
    xlim([rt(1) rt(2)])
    %     ylim([rf(1) rf(2)])
    shading interp
    view(Az,El)
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

%% % graficar areas de ERD en unas frecuencias de un sujeto y un canal % %% Firma
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% grafica las primeras 6 frecuencias del ERD calculado de un sujeto en
% especifico
sub = 3;
ch = 10;
for i = 16
    figure('name',['Sujeto: ' num2str(sub) ' Ch: Cz'],'Tag',['Sujeto: ' num2str(sub) ' freq: ' num2str(f(i))])
    area(t,squeeze(ERD{sub}(ch,i,:)))
    title(['Sujeto: ' num2str(sub) ' Ch: Cz'])
    xlim([1 5])
end

%% %%%%%%%%%%% ERD para mejores frecuencias max en potencia %%%%%%%%%%% %%
% selecciona la mejor deñal segun la maxima potencia calculada en ERD.
ERD_ch = cell(numel(ERD),1);
ff = cell(numel(ERD),1);
for i = 1:numel(ERD)
    % potencia de ERD.
    ERD_tem = sum(ERD{i}(:,5:12,91:end).^2,3);
    % señecciona la mejor potencia de la señal.
    for ch = 1:size(ERD{1},1)
        dd(ch,1) = squeeze(max(ERD_tem(ch,:)));
    end
    for j = 1:size(ERD{1},1)
        ch_tem = find(ERD_tem(j,:) == dd(j));
        if ch_tem < 5
            ch_tem = 5;
        end
        freq_ch(j,1) = ch_tem;
    end
    ff{i} = freq_ch';
    for k = 1:numel(freq_ch)
        ERD_ch{i}(k,:) = ERD{i}(k,freq_ch(k),:);
    end
end
%% area de figuras
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   % posiciones de las graficas.
% rango de amplitud
Lrang = [-2.5 2.5];                                                         % rango de visualizacion en emplitud.
% graficador en 2D todos los canales
freq = 10; % s5
for sub = 5
    i = 1;
    for ch = 10
        figure(sub+30)
        %         subplot(6,7,posi(i))
        %         imagesc(t(16:end)',f(7:33),squeeze(ERD{ch}(7:33,16:end)),[Lrang(1) Lrang(2)])
        area(t(16:end),squeeze(ERD{sub}(ch,freq,16:end)))
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        xlim([1 5])
        %         ylim([1 40]);
        legend(num2str(f(freq)))
        colormap jet
        drawnow
        % if s(sub,1) == 1
        %     title(['Ch: ' num2str(1)])
        % end
        i = i+1;
    end
end
% ff = cell2mat(ff);
%% %%%%%% grafica de ERD para mejores frecuencias max en potencia %%%%% %%   mean and std
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
for sub = 1:9
    i = 1;
    for ch = 1:22
        figure(sub+10)
        subplot(6,7,posi(i))
        channel = squeeze(ERD{sub}(ch,5:20,16:end));
        errorbar(t(16:end),mean(channel),std(channel))
        %         area(t(16:end),squeeze(ERD_ch{sub}(ch,16:end)))
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        ylabel(['Freq' num2str(f(ff{sub}(ch)))])
        xlim([1 4.5])
        drawnow
        i = i+1;
    end
end

%% %%%%%% plot 3 - Grafica de primeras frecuencias de cada sujeto %%%%%% %%
ttt = t; % copia del tiempo para recortar la grafica
% bins time , freq
bins = 5:25;
% Angle graphic.
% Az = -33;
% El = 66;
Az = 2;
El = 1;
ch = 10;  % channel
% plotear
for s =3
    figure(s+20)
    hold on
    for kl = 1:numel(bins)
        ttt(:) = f(bins(kl));
        plot3(t(16:116)',ttt(16:116),squeeze(ERD{s}(ch,kl,16:116)))
    end
    plot3(t.*0+2,t.*6,t.*0,'--rs','LineWidth',2)
    % datos de la grafica
    title(['Subject: ' num2str(s) ' Ch: ' num2str(ch)])
    xlabel('time (seg)')
    ylabel('Freq')
    zlabel('Amp')
    %     xlim([0 7])
    grid on
    view(Az,El);
    hold off
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%% signal TOPOPLOT %%%%%%%%%%%%%%%%%%%%%%%%%% %% --> modificado con movie
% y = ERD_ch; % SEÑAL
% etiquetas 10-20%
load labels
chan = labels;
load layout
[~,pos] = ismember(chan,labels);
% % CSD for the EEG layout
M1.lab = labels(pos);
M1.xy = data(pos,[4 5]);
% tenemos la señal
% control de video.
MovMt = cell(9,1);
for sub = 3 % 1:9
    hold on
    set(0,'DefaultFigureWindowStyle','docked')
    for ufreq = 8
        figure(ufreq-7)
        for utemp = 49
            MyTopo_fun(ERD{sub}(:,ufreq,utemp),M1.xy,M1.lab)
            title(['Sujeto: ' num2str(sub) ' Freq: ' num2str(f(ufreq)) ' Temp: ' num2str(t(utemp))])
            axis square
            shading interp
            colormap jet
            drawnow
            MT(utemp-46) = getframe;
            pause(1)
        end
        MTtemp(ufreq-7,:) = MT;
    end
    hold off
    MovMt{sub} = MTtemp;
end
%%
movie(MovMt{1},26,7)                            % reproduce las imagenes.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% relevance ERD %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
load('relevance.mat')                           % Datos relacionados a la relavancia de los canales.
% relevancia por sujetos y canales
rel = relevance.R(:,1);                         % organiza la relevancia segun los canales
rel = cell2mat(rel);                            % organiza en una sola matriz.
rel = rel./sum(rel,2);                          % relevancia promediada
%% Grafica de solo el mejor canal segun la relevancia. 3D %%
% Plot individual 3D - cada sujeto - mejor canal
% coloca todas las graficas en un mismo conjunto.
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posición de las graficas
posi = [1];
% Rango de freq.
frang = [1 40];
% Rango de tiempo
trang = [16 116];
% rango de graficar limites
rt = [1 5]; % tiempo
rf = [1 30];% freq
% sujetos
% s = [3 2];
% graficador en 3D todos los canales
% sub = 1; % numero de la figura
for sub = 1:9
    % graficas
    figure(sub)
    subplot(1,1,posi(1))
    relmax = max(rel(sub,:));
    ChanCz = find(rel(sub,:) == relmax);
    surf(t(trang(1):trang(2)),f(frang(1):frang(2)),squeeze(ERD{sub}(ChanCz,frang(1):frang(2),trang(1):trang(2))))
    title(['CH: ' labels{ChanCz} ' - Sujeto: ' num2str(sub)])
    %     title(['ch: ' num2str(ChanCz) ' Sujeto: ' num2str(sub)])
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
%% %%%%%%%%%%%%%%%%% Topoplot según la relevancia %%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% signal TOPOPLOT %%%%%%%%%%%%%%%%%%%%%%%%%% %%  --> modificado con movie
ERD_ch = cell(9,1);
for s = 1:9
    for ch = 1:22
        ERD_ch{s}(ch,:,:) = ERD{s}(ch,:,:).*rel(s,ch); % SEÑAL
    end
end
%%
% y = ERD_ch;                                             % SEÑAL
% etiquetas 10-20%
chan = labels;
[~,pos] = ismember(chan,labels);
% % CSD for the EEG layout
M1.lab = labels(pos);
M1.xy = data(pos,[4 5]);
% tenemos la señal
% control de video.
MovMt = cell(9,1);
for sub = 1:9 % 1:9
    hold on
    set(0,'DefaultFigureWindowStyle','docked')
    for ufreq = 7:11
        figure(ufreq-6)
        for utemp = 41:91
            MyTopo_fun(ERD{sub}(:,ufreq,utemp),M1.xy,M1.lab)
            title(['Sujeto: ' num2str(sub) ' Freq: ' num2str(f(ufreq)) ' Temp: ' num2str(t(utemp))])
            axis square
            shading interp
            colormap jet
            drawnow
            MT(utemp-40) = getframe;
            pause(1)
        end
        close all
        MTtemp(ufreq-6,:) = MT;
    end
    hold off
    MovMt{sub} = MTtemp;
end
%%
movie(MovMt{1},26,7)                                                        % reproduce las imagenes.
%% Selecciona el mejor canal segun las mejores señales ERD
ff = cell2mat(ff);
ERDchm = cell(9,1);
for s =1:9
    for ch = 1:22
        dat = ff(s,ch);
        ERDchm{s}(ch,:) = squeeze(ERD_ch{s}(dat,16:end)).*rel(s,ch);
    end
end

%% ERD para mejores frecuencias mejor con entropia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ERD_ch = cell(numel(ERD),1);
ff = cell(numel(ERD),1);
for i = 1:numel(ERD)
    ERD_tem = sum(ERD{i}(:,5:17,91:end).^2,3);
    for ch = 1:size(ERD{1},1)
        dd(ch,1) = squeeze(max(ERD_tem(ch,:)));
    end
    for j = 1:size(ERD{1},1)
        ch_tem = find(ERD_tem(j,:) == dd(j));
        freq_ch(j,1) = ch_tem;
    end
    ff{i} = freq_ch;
    for k = 1:numel(freq_ch)
        ERD_ch{i}(k,:) = ERD{i}(k,freq_ch(k),:);
    end
end

%% %%%% grafica de ERD para mejores frecuencias mejor con entropia %%%% %% -----> para probar con un proceso diferente
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% posiciones de las graficas
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];   %
% ciclo de graficas.
for sub = 1:9
    i = 1;
    for ch = 1:22
        figure(sub+10)
        subplot(6,7,posi(i))
        area(t(16:end),squeeze(ERD_ch{sub}(ch,16:end)))
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        ylabel(['Freq' num2str(f(ff{sub}(ch)))])
        xlim([0 7])
        drawnow
        i = i+1;
    end
end

%% plot 3
ttt = t;
% bins time , freq
bins = 4:33;
% Angle graphic.
Az = -33;
El = 66;
ch = 10;                                                                    % channel seleccionado
for s =1:9
    figure(s+20)                    %
    hold on
    for kl = 1:numel(bins)
        ttt(:) = f(bins(kl));
        plot3(t(16:end)',ttt(16:end),squeeze(ERD{s}(ch,kl,16:end)))
    end
    plot3(t.*0+2,t.*6,t.*0,'--rs','LineWidth',2)
    title(['Subject: ' num2str(s) ' Ch: ' num2str(ch)])
    xlabel('time (seg)')
    ylabel('Freq')
    zlabel('Amp')
    xlim([0 7])
    grid on
    view(Az,El);
    hold off
end

%% Headplot for channels

% load ch_sel
%
for sub = 1:9
    %     figure(sub)
    %     subplot(6,7,4)
    %     imagesc(t',f',squeeze(ERD{sub}(1,:,:)),[-1 1])
    %     if s(sub,1) == 1
    %     title(['Ch: ' num2str(1)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %     end
    %
    %     subplot(6,7,9)
    %     imagesc(t',f',squeeze(ERD{sub}(2,:,:)),[-1 1])
    %     if s(sub,2) == 1
    %     title(['Ch: ' num2str(2)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,10)
    %     imagesc(t',f',squeeze(ERD{sub}(3,:,:)),[-1 1])
    %     if s(sub,3) == 1
    %     title(['Ch: ' num2str(3)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,11)
    %     imagesc(t',f',squeeze(ERD{sub}(4,:,:)),[-1 1])
    %     if s(sub,4) == 1
    %     title(['Ch: ' num2str(4)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,12)
    %     imagesc(t',f',squeeze(ERD{sub}(5,:,:)),[-1 1])
    %     if s(sub,5) == 1
    %     title(['Ch: ' num2str(5)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,13)
    %     imagesc(t',f',squeeze(ERD{sub}(6,:,:)),[-1 1])
    %     if s(sub,6) == 1
    %     title(['Ch: ' num2str(6)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,15)
    %     imagesc(t',f',squeeze(ERD{sub}(7,:,:)),[-1 1])
    %     if s(sub,7) == 1
    %     title(['Ch: ' num2str(7)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,16)
    %     imagesc(t',f',squeeze(ERD{sub}(8,:,:)),[-1 1])
    %     if s(sub,8) == 1
    %     title(['Ch: ' num2str(8)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,17)
    %     imagesc(t',f',squeeze(ERD{sub}(9,:,:)),[-1 1])
    %     if s(sub,9) == 1
    %     title(['Ch: ' num2str(9)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,18)
    %     imagesc(t',f',squeeze(ERD{sub}(10,:,:)),[-1 1])
    %     if s(sub,10) == 1
    %     title(['Ch: ' num2str(10)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,19)
    %     imagesc(t',f',squeeze(ERD{sub}(11,:,:)),[-1 1])
    %     if s(sub,11) == 1
    %     title(['Ch: ' num2str(11)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,20)
    %     imagesc(t',f',squeeze(ERD{sub}(12,:,:)),[-1 1])
    %     if s(sub,12) == 1
    %     title(['Ch: ' num2str(12)])
    %     xlabel('time (seg)')
    %
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,21)
    %     imagesc(t',f',squeeze(ERD{sub}(13,:,:)),[-1 1])
    %     if s(sub,13) == 1
    %     title(['Ch: ' num2str(13)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,23)
    %     imagesc(t',f',squeeze(ERD{sub}(14,:,:)),[-1 1])
    %     if s(sub,14) == 1
    %     title(['Ch: ' num2str(14)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,24)
    %     imagesc(t',f',squeeze(ERD{sub}(15,:,:)),[-1 1])
    %     if s(sub,15) == 1
    %     title(['Ch: ' num2str(15)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,25)
    %     imagesc(t',f',squeeze(ERD{sub}(16,:,:)),[-1 1])
    %     if s(sub,16) == 1
    %     title(['Ch: ' num2str(16)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,26)
    %     imagesc(t',f',squeeze(ERD{sub}(17,:,:)),[-1 1])
    %     if s(sub,17) == 1
    %     title(['Ch: ' num2str(17)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,27)
    %     imagesc(t',f',squeeze(ERD{sub}(18,:,:)),[-1 1])
    %     if s(sub,18) == 1
    %     title(['Ch: ' num2str(18)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,31)
    %     imagesc(t',f',squeeze(ERD{sub}(19,:,:)),[-1 1])
    %     if s(sub,19) == 1
    %     title(['Ch: ' num2str(19)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,32)
    %     imagesc(t',f',squeeze(ERD{sub}(20,:,:)),[-1 1])
    %     if s(sub,20) == 1
    %     title(['Ch: ' num2str(20)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,33)
    %     imagesc(t',f',squeeze(ERD{sub}(21,:,:)),[-1 1])
    %     if s(sub,21) == 1
    %     title(['Ch: ' num2str(21)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    %     subplot(6,7,39)
    %     imagesc(t',f',squeeze(ERD{sub}(22,:,:)),[-1 1])
    %     if s(sub,22) == 1
    %     title(['Ch: ' num2str(22)])
    %     xlabel('time (seg)')
    %     ylim([1 40]);
    %     colormap jet
    %     drawnow
    %
    %     end
    %
    % %     label_save = ['Canal_' num2str(sub)];
    % %     saveas(eval('fig1'),[label_save '.fig'])
end

%% Plot ERD all subjects, all channels
r=[1 2 3];
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
for sub = 3:3
    for ch = 1:22
        figure(ch)
        imagesc(t',f',squeeze(ERD{sub}(ch,:,:)),[-1 1])
        title(['Ch: ' num2str(ch)])
        ylabel('Frequency (Hz)')
        xlabel('time (seg)')
        ylim([1 40]);
        colormap jet
        drawnow
    end
    pause(10)
    c = 22;
    figure(23)
    plot(t,squeeze(ERD{sub}(c,8,:)))
    figure(24)
    plot(t,squeeze(ERD{sub}(c,8,:)))
end
% figure(25)
% plot(t,squeeze(ERD{sub}(c,9,:)))

%% %%%%%%%%%%%%%%%%%%%%%%%%% Plot spectrogram  %%%%%%%%%%%%%%%%%%%%%%%%% %% -
for ch = 1:22
    figure(ch)
    clf
    P = abs(log(squeeze(Sujetos{8}(50,ch,:,:))));
    %imagesc(time',f',P)
    contourf(t',f',P,30)
    colormap parula
    title(['Sujeto: ' num2str(sub)])
    ylabel('Frequency (Hz)')
    xlabel('time (seg)')
    ylim([1 40]);
    set(gca,'YScale','log')
    drawnow
end

%% Seccion B
%% %%%%%%%%%%%%%%%%%%%%%%% ERDs - WAVELET PACKET %%%%%%%%%%%%%%%%%%%%%%% %%
%% Database
clear all; clc % limpiar datos.

% load data
load('BCICIV_2a.mat');

% ERD con Wavelet band
t = 0:1/250:7-(1/250);
t1 = 0;
t2 = 2;
lv = 5;
subjeto = 3;
name = 'ad5';
for su = subjeto%1:9
    Sujeto = X{su};
    etiqueta = y{su};
    
    % escojer los train de las etiquetas 1 y 2
    Data = Sujeto(ismember(etiqueta,[2]));
    
    % escojer el que tenga la menor cantidad de trials
    N_trial = length(Data);
    Canales = [1:22]; N_canal = length(Canales);
    
    for tri = 1:N_trial
        fprintf(['sujeto:1 ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        % signal
        signal = Data{tri}(:,:);
        % function filterwavelet
        [X_dat, T] = FilterWave(signal,'sym2',lv);
        Y_suj{su}(tri,:,:,:) = X_dat;
    end
end
% tiempo de referencia
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - t2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
for su = subjeto
    fprintf(['sujeto: ' num2str(su) ' de ' '9' '\n'])
    % Et{Pnc(t,f):te[Ta,Tb]}
    r_nc{su,1} = squeeze(mean(Y_suj{su}(:,:,:,ul:up),4));
end
fprintf(['En{r_nc(f)}\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{r_nc(f)}
    r_c{y,1} = squeeze(mean(r_nc{y,1},1));
end
fprintf(['En{Pnc(t,f)}\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % En{Pnc(t,f)}
    m_c{y,1} = squeeze(mean(Y_suj{y},1));
end
fprintf(['ERD = m_c/r_c\n'])
for y = subjeto
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1;
end

% Headplot for channels (lv = 1:n)

sub = subjeto;
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
if lv == 4
    hh = [1 3 7 15];                                                        % 31 %34 8 17 35 36 18 37 38 4 9 19 39 40 20 41 42 10 21 43 44 % hh = 1:62;
elseif lv == 5
    hh = [1 3 7 15 31];
elseif lv == 6
    hh = [1 3 7 15 31 63];
end
f = 1:length(hh);
da = zeros(length(hh),22,1750);

for i =1:length(hh)
    da(i,:,:) = ERD{sub}(hh(i),:,:);
end
% imagesc(t',f,squeeze(da(:,ch,:)),[l1 l2])
% Headplot for all_channels
l1 = -5;
l2 = 5;
% -2.2359e+03 2.2359e+03
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
i = 1;
for ch = 1:22
    figure(subjeto)
    subplot(6,7,posi(i))
    imagesc(t',f,squeeze(da(:,ch,:)),[l1 l2])
    title(['Ch: ' num2str(ch)])
    i = i+1;
end

%% Headplot for signals in one channel
set(0,'DefaultFigureWindowStyle','docked')                                  % grafica cada figura en una misma pantalla.
% for j = 7:13
for i = 1:length(hh)
    figure(i+10)
    plot((1:1750)/250,squeeze(da(i,12,:)))
end
pause
% end
%% Save (data)
save(name,'f','t','da','hh')
%% Headplot for 1_channel
% % for ch = 1:22
%     figure(22)
% %         subplot(6,7,4)
% %     subplot(6,7,4)
%     imagesc(t',f,squeeze(da(:,22,:)),[l1 l2])
%     title(['Ch: ' num2str(22)])
% end
%% plot - subplot individual
%         subplot(6,7,4)
%         imagesc(t',f,squeeze(da(:,1,:)),[l1 l2])
%         subplot(6,7,9)
%         imagesc(t',f,squeeze(da(:,2,:)),[l1 l2])
% %         if s(sub,2) == 1
% %         title(['Ch: ' num2str(2)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,10)
%         imagesc(t',f,squeeze(da(:,3,:)),[l1 l2])
% %         if s(sub,3) == 1
% %         title(['Ch: ' num2str(3)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,11)
%         imagesc(t',f,squeeze(da(:,4,:)),[l1 l2])
% %         if s(sub,4) == 1
% %         title(['Ch: ' num2str(4)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,12)
%         imagesc(t',f,squeeze(da(:,5,:)),[l1 l2])
% %         if s(sub,5) == 1
% %         title(['Ch: ' num2str(5)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,13)
%         imagesc(t',f,squeeze(da(:,6,:)),[l1 l2])
% %         if s(sub,6) == 1
% %         title(['Ch: ' num2str(6)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,15)
%         imagesc(t',f,squeeze(da(:,7,:)),[l1 l2])
% %         if s(sub,7) == 1
% %         title(['Ch: ' num2str(7)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,16)
%         imagesc(t',f,squeeze(da(:,8,:)),[l1 l2])
% %         if s(sub,8) == 1
% %         title(['Ch: ' num2str(8)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,17)
%         imagesc(t',f,squeeze(da(:,9,:)),[l1 l2])
% %         if s(sub,9) == 1
% %         title(['Ch: ' num2str(9)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,18)
%         imagesc(t',f,squeeze(da(:,10,:)),[l1 l2])
% %         if s(sub,10) == 1
% %         title(['Ch: ' num2str(10)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,19)
%         imagesc(t',f,squeeze(da(:,11,:)),[l1 l2])
% %         if s(sub,11) == 1
% %         title(['Ch: ' num2str(11)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,20)
%         imagesc(t',f,squeeze(da(:,12,:)),[l1 l2])
% %         if s(sub,12) == 1
% %         title(['Ch: ' num2str(12)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,21)
%         imagesc(t',f,squeeze(da(:,13,:)),[l1 l2])
% %         if s(sub,13) == 1
% %         title(['Ch: ' num2str(13)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,23)
%         imagesc(t',f,squeeze(da(:,14,:)),[l1 l2])
% %         if s(sub,14) == 1
% %         title(['Ch: ' num2str(14)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,24)
%         imagesc(t',f,squeeze(da(:,15,:)),[l1 l2])
% %         if s(sub,15) == 1
% %         title(['Ch: ' num2str(15)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,25)
%         imagesc(t',f,squeeze(da(:,16,:)),[l1 l2])
% %         if s(sub,16) == 1
% %         title(['Ch: ' num2str(16)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,26)
%         imagesc(t',f,squeeze(da(:,17,:)),[l1 l2])
% %         if s(sub,17) == 1
% %         title(['Ch: ' num2str(17)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,27)
%         imagesc(t',f,squeeze(da(:,18,:)),[l1 l2])
% %         if s(sub,18) == 1
% %         title(['Ch: ' num2str(18)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,31)
%         imagesc(t',f,squeeze(da(:,19,:)),[l1 l2])
% %         if s(sub,19) == 1
% %         title(['Ch: ' num2str(19)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,32)
%         imagesc(t',f,squeeze(da(:,20,:)),[l1 l2])
% %         if s(sub,20) == 1
% %         title(['Ch: ' num2str(20)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%         subplot(6,7,33)
%         imagesc(t',f,squeeze(da(:,21,:)),[l1 l2])
% %         if s(sub,21) == 1
% %         title(['Ch: ' num2str(21)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
%         subplot(6,7,39)
%         imagesc(t',f,squeeze(da(:,22,:)),[l1 l2])
% %         if s(sub,22) == 1
% %         title(['Ch: ' num2str(22)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %         end
%     label_save = ['Canal_' num2str(sub)];
%     saveas(eval('fig1'),[label_save '.fig'])
%% prueba plot
% for i = 1:22
% cnl = i;
%
% f_plot(Sujetos,cnl,paso,Nw)
%
%     canal = i; paso = Nw-noverlap;
%
%     label_save = ['Canal_' num2str(canal)];
%     saveas(eval('fig1'),[label_save '.png'])
%     saveas(eval('fig1'),[label_save '.fig'])
% end
%% plot en coordenada
% fig = gca;
% subplot(fig,'Position',[0,0,0.5,1])
% plot(ERD{1}(1,5,:))