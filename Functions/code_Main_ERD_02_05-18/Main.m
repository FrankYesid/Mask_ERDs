%% Programa ERD Spectrogram & Wavelet
% 
% F. Y. Zapata C.
%% Seccion A
%% %%%%%%%%%%%%%%%%%%%%%%%% ERDs - SPECTROGRAM %%%%%%%%%%%%%%%%%%%%%%%% %%
%% Database 
% clear all; clc
%% LOAD DATA
clear all;clc
load('BCICIV_2a.mat');
% ERDs spectrogram
% parametros de la STFT
% Nw = 1; % tama?o de la ventana en seg
Nw = 0.8*fs;
window = hamming(Nw);
noverlap = floor(0.95*Nw);

X_suj = cell(9,1);

for su = 1:9
    
    Sujeto = X{su};
    etiqueta = y{su};
    
    % escojer los train de las etiquetas 1 y 2
    Izq_Der = Sujeto(ismember(etiqueta,[1 2]));
    
    % escojer el que tenga la menor cantidad de trials
    N_trial = length(Izq_Der);
    Canales = [1:22]; N_canal = length(Canales);
    
    % Inicializar matrices
    temp = spectrogram(Izq_Der{1}(:,1), window, noverlap,Nw);
    temp = size(temp,2);
 
    X_suj(su) = {zeros(N_trial,N_canal,floor((Nw/2)+1),temp)}; % Trials x canal x frecuencia x tiempo
    
    for tri = 1:N_trial
        fprintf(['sujeto: ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        for cnl = 1:N_canal
            % signal
            signal = Izq_Der{tri}(:,Canales(cnl));
            
            % Calcular STFT
            [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
            X_Class = abs(X_Class);
            
            % Almacenar espectrogramas
            X_suj{su}(tri,cnl,:,:) = X_Class;            
        end
    end
end

Sujetos = X_suj;
t1 = 0;
t2 = 2;
temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - 2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
r_nc = cell(9,1);
r_c = cell(9,1);
m_c = cell(9,1);
fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])

for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % Et{Pnc(t,f):te[Ta,Tb]}
    r_nc{y,1} = squeeze(mean(X_suj{y}(:,:,:,ul:up),4));
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
fprintf(['ERD = m_c/r_c\n'])
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1;
end

% % c = 22;
% % figure(23)
% % plot(t,squeeze(ERD{sub}(c,8,:)))
% % save('ERD_t13.mat','ERD')

%%  Analisis de relevancia -------------NO---------------------------------
% % load('X_sujetos_Izq_Der.mat')
% % 
% % Sujetos = f_pca(X_sujetos_class_1);
% % 
% % save('Sujetos.mat','Sujetos');

%% Graficar ---------------------------------------------------------------
% 
% % clc
% % fs = 250;
% % Nw = 1; % tama?o de la ventana en seg
% % Nw = 1*fs;
% % window = hamming(Nw);
% % noverlap = floor(0.9*fs);
% 
% % load('ERD.mat')
% % paso = Nw-noverlap;
% 
% % load('Sujetos.mat')
% % [N_freq N_t] = size(squeeze(Sujetos{1}(1,1,:,:)));
% % f = [0:floor((N_freq-1))]*fs/(2*N_freq);
% % time = [0:N_t-1]*paso/fs+(Nw/2)/fs;

%% Headplot for channel ciclo
l1 = -1;
l2 = 1;
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
for sub = 1:9
    i = 1;
    for ch = 1:22
        figure(sub)
        subplot(6,7,posi(i))
        imagesc(t',f',squeeze(ERD{sub}(ch,:,:)),[l1 l2])
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        xlim([0 7])
        ylim([1 40]);
        colormap jet
        drawnow
        % if s(sub,1) == 1
        %     title(['Ch: ' num2str(1)])
        % end
        i = i+1;
    end
end

%% graficar areas de ERD en unas frecuencias
set(0,'DefaultFigureWindowStyle','docked')
for i = 1:6
    figure(20)
    subplot(3,2,i)
    area(t,squeeze(ERD{2}(20,i,:)))
end

%% grafica de ERD para mejores frecuencias
ERD_ch = cell(numel(ERD),1);
ff = cell(numel(ERD),1);
for i = 1:numel(ERD)
    ERD_tem = sum(ERD{i}(:,:,1:8).^2,3);
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
    ff{i} = freq_ch;
    for k = 1:numel(freq_ch)
        ERD_ch{i}(k,:) = ERD{i}(k,freq_ch(k),:);
    end
end
%%
posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
for sub = 1:9
    i = 1;
    for ch = 1:22
        figure(sub+10)
        subplot(6,7,posi(i))
        area(t,squeeze(ERD_ch{sub}(ch,:)))
        title(['Ch: ' num2str(ch)])
        xlabel('time (seg)')
        ylabel(['Freq' num2str(f(ff{sub}(ch)))])
        xlim([0 7])
        drawnow
        % if s(sub,1) == 1
        %     title(['Ch: ' num2str(1)])
        % end
        i = i+1;
    end
end

%% Headplot for channels

% load ch_sel
%
% for sub = 1:9
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
% end

%% Plot ERD all subjects, all channels
% % r=[1 2 3];
% % set(0,'DefaultFigureWindowStyle','docked')
% % for sub = 3:3
% %     for ch = 1:22
% %         figure(ch)
% %         imagesc(t',f',squeeze(ERD{sub}(ch,:,:)),[-1 1])
% %         title(['Ch: ' num2str(ch)])
% %         ylabel('Frequency (Hz)')
% %         xlabel('time (seg)')
% %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% %          
% %     end
% %     pause(10)
% % 
% % c = 22;
% % figure(23)
% % plot(t,squeeze(ERD{sub}(c,8,:)))
% % 
% % figure(24)
% % plot(t,squeeze(ERD{sub}(c,8,:)))
% % 
% % end
% % % figure(25)
% % % plot(t,squeeze(ERD{sub}(c,9,:)))

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
% %% Database 
% clear all; clc % limpiar datos.
% 
% % load data
% load('BCICIV_2a.mat');
% 
% % ERD con Wavelet band
% t = 0:1/250:7-(1/250);
% t1 = 0;
% t2 = 2;
% lv = 5;
% subjeto = 9;
% name = 'ad5';
% for su = subjeto%1:9
%     Sujeto = X{su};
%     etiqueta = y{su};
%     
%     % escojer los train de las etiquetas 1 y 2
%     Data = Sujeto(ismember(etiqueta,[2]));
%     
%     % escojer el que tenga la menor cantidad de trials
%     N_trial = length(Data);
%     Canales = [1:22]; N_canal = length(Canales);
%     
%     for tri = 1:N_trial
%         fprintf(['sujeto:1 ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
%         % signal
%         signal = Data{tri}(:,:);
%         % function filterwavelet 
%         [X_dat, T] = FilterWave(signal,'sym2',lv);
%         Y_suj{su}(tri,:,:,:) = X_dat; 
%     end     
% end
% 
% 
% % tiempo de referencia
% temp1 = abs(t - t1);
% min1 = min(temp1);
% temp2 = abs(t - t2);
% min2 = min(temp2);
% ul = find(temp1 == min1);
% up = find(temp2 == min2);
% 
% %
% fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
% for su = subjeto
%     fprintf(['sujeto: ' num2str(su) ' de ' '9' '\n'])
%     % Et{Pnc(t,f):te[Ta,Tb]}
%     r_nc{su,1} = squeeze(mean(Y_suj{su}(:,:,:,ul:up),4));
% end
% 
% %
% fprintf(['En{r_nc(f)}\n'])
% for y = subjeto
%     fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
%     % En{r_nc(f)}
%     r_c{y,1} = squeeze(mean(r_nc{y,1},1));
% end
% %
% fprintf(['En{Pnc(t,f)}\n'])
% for y = subjeto
%     fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
%     % En{Pnc(t,f)}
%     m_c{y,1} = squeeze(mean(Y_suj{y},1));
% end
% 
% fprintf(['ERD = m_c/r_c\n'])
% for y = subjeto
%     fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
%     % ERD = m_c/r_c
%     ERD{y,1}= bsxfun(@times,m_c{y,1},1./r_c{y,1}) - 1;
% end
% 
% % Headplot for channels (lv = 1:n)
% 
% sub = subjeto;
% set(0,'DefaultFigureWindowStyle','docked')
% if lv == 4
%     hh = [1 3 7 15]; % 31 %34 8 17 35 36 18 37 38 4 9 19 39 40 20 41 42 10 21 43 44 % hh = 1:62;
% elseif lv == 5
%     hh = [1 3 7 15 31];
% elseif lv == 6
%     hh = [1 3 7 15 31 63];
% end
%  
% f = 1:length(hh);
% da = zeros(length(hh),22,1750);
% 
% for i =1:length(hh)
%         da(i,:,:) = ERD{sub}(hh(i),:,:);
% end
% 
% % imagesc(t',f,squeeze(da(:,ch,:)),[l1 l2])
% 
% % Headplot for all_channels
% l1 = -5;
% l2 = 5;
% % -2.2359e+03 2.2359e+03
% posi = [4 9 10 11 12 13 15 16 17 18 19 20 21 23 24 25 26 27 31 32 33 39];
% i = 1;
% for ch = 1:22
%     figure(subjeto)
%     subplot(6,7,posi(i))
%     imagesc(t',f,squeeze(da(:,ch,:)),[l1 l2])
%     title(['Ch: ' num2str(ch)])
%     i = i+1;
% end
% 
% %% Headplot for signals in one channel
% set(0,'DefaultFigureWindowStyle','docked')
% % for j = 7:13
%     for i = 1:length(hh)
%         figure(i+10)
%         plot((1:1750)/250,squeeze(da(i,12,:)))
%     end
%    pause
% % end
% %% Save (data) 
% save(name,'f','t','da','hh')

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
% % 
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
% 
% %         end
% 
%         subplot(6,7,10)
%         imagesc(t',f,squeeze(da(:,3,:)),[l1 l2])
% %         if s(sub,3) == 1
% %         title(['Ch: ' num2str(3)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,11)
%         imagesc(t',f,squeeze(da(:,4,:)),[l1 l2])
% %         if s(sub,4) == 1
% %         title(['Ch: ' num2str(4)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,12)
%         imagesc(t',f,squeeze(da(:,5,:)),[l1 l2])
% %         if s(sub,5) == 1
% %         title(['Ch: ' num2str(5)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,13)
%         imagesc(t',f,squeeze(da(:,6,:)),[l1 l2])
% %         if s(sub,6) == 1
% %         title(['Ch: ' num2str(6)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,15)
%         imagesc(t',f,squeeze(da(:,7,:)),[l1 l2])
% %         if s(sub,7) == 1
% %         title(['Ch: ' num2str(7)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,16)
%         imagesc(t',f,squeeze(da(:,8,:)),[l1 l2])
% %         if s(sub,8) == 1
% %         title(['Ch: ' num2str(8)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,17)
%         imagesc(t',f,squeeze(da(:,9,:)),[l1 l2])
% %         if s(sub,9) == 1
% %         title(['Ch: ' num2str(9)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,18)
%         imagesc(t',f,squeeze(da(:,10,:)),[l1 l2])
% %         if s(sub,10) == 1
% %         title(['Ch: ' num2str(10)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,19)
%         imagesc(t',f,squeeze(da(:,11,:)),[l1 l2])
% %         if s(sub,11) == 1
% %         title(['Ch: ' num2str(11)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,20)
%         imagesc(t',f,squeeze(da(:,12,:)),[l1 l2])
% %         if s(sub,12) == 1
% %         title(['Ch: ' num2str(12)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,21)
%         imagesc(t',f,squeeze(da(:,13,:)),[l1 l2])
% %         if s(sub,13) == 1
% %         title(['Ch: ' num2str(13)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,23)
%         imagesc(t',f,squeeze(da(:,14,:)),[l1 l2])
% %         if s(sub,14) == 1
% %         title(['Ch: ' num2str(14)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,24)
%         imagesc(t',f,squeeze(da(:,15,:)),[l1 l2])
% %         if s(sub,15) == 1
% %         title(['Ch: ' num2str(15)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,25)
%         imagesc(t',f,squeeze(da(:,16,:)),[l1 l2])
% %         if s(sub,16) == 1
% %         title(['Ch: ' num2str(16)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,26)
%         imagesc(t',f,squeeze(da(:,17,:)),[l1 l2])
% %         if s(sub,17) == 1
% %         title(['Ch: ' num2str(17)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,27)
%         imagesc(t',f,squeeze(da(:,18,:)),[l1 l2])
% %         if s(sub,18) == 1
% %         title(['Ch: ' num2str(18)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,31)
%         imagesc(t',f,squeeze(da(:,19,:)),[l1 l2])
% %         if s(sub,19) == 1
% %         title(['Ch: ' num2str(19)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,32)
%         imagesc(t',f,squeeze(da(:,20,:)),[l1 l2])
% %         if s(sub,20) == 1
% %         title(['Ch: ' num2str(20)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
%         subplot(6,7,33)
%         imagesc(t',f,squeeze(da(:,21,:)),[l1 l2])
% %         if s(sub,21) == 1
% %         title(['Ch: ' num2str(21)])
% %         xlabel('time (seg)')
% %         xlim([0 7])
% % %         ylim([1 40]);
% %         colormap jet
% %         drawnow
% % 
% %         end
% 
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
% 


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

%
%% plot en coordenada
% fig = gca;
% subplot('Position',[0,0,0.5,1])
% plot(ERD{1}(1,5,:))
