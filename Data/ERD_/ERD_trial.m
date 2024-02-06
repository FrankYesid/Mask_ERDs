function [ERD_tr] = ERD_trial(X,freq,tref,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  The function calculates the ERDs of data series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% X: datos del tiempo de la señal. (canales x tiempo)
% freq: rango de frecuencias.
% tref: rango de tiempos de referencia. 
% T: tiempo para analisis para el mejor erd.
% -------------------------------------------------------------------------
% Output:
% ERD_tr: datos del tiempo de la señal. (canales x tiempo). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram.
Nw = 0.8*fs;                                                                % tamaño de la ventana.
window = hamming(Nw);                                                       % Vector de la ventana.
noverlap = floor(0.95*Nw);                                                  % interpolación de la ventana.
% calculo del espectrograma relacionada a EEG.
Izq_Der = X;                                                                % seleccion segun la clase o clases seleccionadas.
Canales = [1:22]; N_canal = length(Canales);                                % cantidad de canales.
% Inicializar matrices
temp = spectrogram(Izq_Der(:,1), window, noverlap,Nw);                      % tamaño temporal de los datos en cada sujeto.
temp = size(temp,2);
X_suj = zeros(N_canal,floor((Nw/2)+1),temp);                                % canal x frecuencia x tiempo.
for cnl = 1:N_canal
    % signal
    signal = Izq_Der(:,Canales(cnl));                                       % Señal de cada canal.
    % Calcular STFT
    [X_Class, f, t] = spectrogram(signal, window, noverlap,Nw,fs);          % Espectrograma de la señal.
    X_Class = abs(X_Class);                                                 % absoluto para tener la señal (freq,tiempo)
    % Almacenar espectrogramas
    X_suj(cnl,:,:) = X_Class;                                               % canal x frecuencia x tiempo
end
% tiempo de referencia para el calculo del ERD.
temp1 = abs(t - tref(1));        
min1 = min(temp1);
temp2 = abs(t - tref(2));
min2 = min(temp2);
ul = find(temp1 == min1);
up = find(temp2 == min2);
r_c = zeros(size(X_suj,1),size(X_suj,2));                                   % celda para esperanza respecto al tiempo de referencia.
ERD_temp = zeros(size(X_suj,1),size(X_suj,2),size(X_suj,3));   
for y = 1:22
    r_c(y,:) = squeeze(mean(X_suj(y,:,ul:up),3));                           % calculo de la referencia.
end
for y = 1:22
    ERD_temp(y,:,:) = bsxfun(@times,X_suj(y,:,:),1./r_c(y,1)) - 1;          % channel,freq,time
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localiza el segmento solicitado.
% frecuencia.
fr1 = abs(f - freq(1));        
min1 = min(fr1);
fr2 = abs(f - freq(2));
min2 = min(fr2);
f1 = find(fr1 == min1);
f2 = find(fr2 == min2);
% tiempo.
temp1 = abs(t - T(1));        
min1 = min(temp1);
temp2 = abs(t - T(2));
min2 = min(temp2);
t1 = find(temp1 == min1);
t2 = find(temp2 == min2);
ERD_ch = cell(22,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ciclos de la señal segun el segmento.
for y = 1:22
    ET = sum(ERD_temp(f1:f2,t1:t2).^2,2);                                   % selecciona freq(8-13) - tiempos(1-5) generalmente.
    % selecciona la mejor potencia de la señal.
    dd = max(ET);
    dat = find(ET == dd) + 6;
    ERD_ch{y} = squeeze(ERD_temp(y,dat,:))';
end
ERD_tr = cell2mat(ERD_ch);