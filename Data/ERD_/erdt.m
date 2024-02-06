clear all; clc
load('BCICIV_2a_All.mat');

%%  data
pos= find(mrkType{1}==768);
%%

% parametros de la STFT
% Nw = 1; % tama?o de la ventana en seg
Nw = 0.7*fs;
window = hamming(Nw);
noverlap = floor(0.9*Nw);

X_suj = cell(9,1);

for su = 1:9
    
    Sujeto = X{su};
    etiqueta = y{su};
    
    % escojer los train de las etiquetas 1 y 2
%     Izq_Der = Sujeto(ismember(etiqueta,[1]));
    
    % escojer el que tenga la menor cantidad de trials
%     N_trial = length(Sujeto);
    Canales = [1:22]; N_canal = length(Canales);
    
    % Inicializar matrices
    temp = spectrogram(Sujeto{1}(:,1), window, noverlap,Nw);
    temp = size(temp,2);
 
    X_suj(su) = {zeros(N_canal,floor((Nw/2)+1),temp)}; % Trials x canal x frecuencia x tiempo
    
%     for tri = 1:N_trial
%         fprintf(['sujeto: ' num2str(su) ' de ' '9' ' ...trial: ' num2str(tri) ' de ' num2str(N_trial) '\n'])
        for cnl = 1:N_canal
            % signal
            signal = Sujeto{su}(:,Canales(cnl));
            
            % Calcular STFT
            [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
            X_Class = abs(X_Class);
            
            % Almacenar espectrogramas
            X_suj{su}(cnl,:,:) = X_Class;            
        end
%     end
end
%%

Sujetos = X_suj;
t1 = 0:7:70;
t2 = 2:7:72;

temp1 = abs(t - t1);
min1 = min(temp1);
temp2 = abs(t - 2);
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);

fprintf(['Et{Pnc(t,f):te[Ta,Tb]}\n'])
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % Et{Pnc(t,f):te[Ta,Tb]}
    r_nc{y,1} = squeeze(mean(X_suj{y}(:,:,:,ul:up),4));
end

% fprintf(['En{r_nc(f)}\n'])
% for y = 1:9
%     fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
%     % En{r_nc(f)}
%     r_c{y,1} = squeeze(mean(r_nc{y,1},1));
% end
% 
% fprintf(['En{Pnc(t,f)}\n'])
% for y = 1:9
%     fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
%     % En{Pnc(t,f)}
%     m_c{y,1} = squeeze(mean(X_suj{y},1));
% end

fprintf(['ERD = m_c/r_c\n'])
for y = 1:9
    fprintf(['sujeto: ' num2str(y) ' de ' '9' '\n'])
    % ERD = m_c/r_c
    ERD{y,1}= bsxfun(@times,X_suj{y,1},1./r_nc{y,1}) - 1;
end

