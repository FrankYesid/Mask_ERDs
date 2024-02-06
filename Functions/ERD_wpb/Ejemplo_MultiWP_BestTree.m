clear all; clc;
%% Multi-signal WPD decompostion
%number of samples
ns=20000;
%number of channels
nc=20;
%multi-channel signal
y=rand(ns,nc);
% load('y')

%% Multi-channel Wavelet Packet
Wave = 'sym15';
lv = 4;
[WaveRec  WaveCoef]= MultiWPdec(X{1}{1},'sym15',4);

%% BestTree: selecionar los nodos mas representativos
% BestRec: contiene los mejores nodos para la reconstruccion de la se?al Waverec(time,fil,canal)
% BestT: Contiene los nodos mas representativos pot nivel BestT{lv}(time,fil,canal) 
[BestRec,BestT] = f_BestTree(WaveRec,WaveCoef);

%% Graficar
y=X{1}{1};
cnl = 1;
A = squeeze(BestRec(:,:,cnl));
y_r = sum(A,2);

figure(1)
Error = sum((y_r-y(:,cnl)).^2)/length(y);
plot(y_r); hold on; plot(y(:,cnl)); hold off; title(['Error: ' num2str(Error)])
legend('Original Signal','Reconstructed signal')

figure(2)
title('Espectro de cada filtro')
fft_plot(A,2)

