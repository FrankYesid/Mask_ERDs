function f_plot(Sujetos,cnl,paso,Nw)
fs = 250;
[N_freq N_t] = size(squeeze(Sujetos{1}(1,:,:)));
f = [0:floor((N_freq-1)/2)]*fs/N_freq;
% time = [0:N_t]*paso/Fs+(Nw/2)/fs ;
time = [0:N_t-1]*paso/fs+(Nw/2)/fs;
for su = 1:9
    subplot(3,3,su)
    imagesc(time',f',Sujetos{su}(cnl,:,:))
    title(['Sujeto: ' num2str(su)])
    ylabel('Frequency (Hz)')
    xlabel('time (seg)')
    ylim([0 40]);
end
suptitle(['Canal ' num2str(cnl)])