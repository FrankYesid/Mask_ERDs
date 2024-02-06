function fft_plot(X,Fs)


L = size(X,1);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:round(L/2)+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
f = Fs*(0:round(L/2))/L;
plot(f,P1) 
xlabel('f (Hz)')
ylabel('|P1(f)|')