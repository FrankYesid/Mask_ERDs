%% Main para graficar la entropias

% entropy
% entropy = rand(2,10,22);
%% parametros
% limites de la amplitud de la señal
limi = [0,1];
% labels del eje X
labelx = num2str([0,1]);
% labels del eje Y
labely = num2str([0,1]);
sel = 1:22;
% limits_y = [0 1];
% limits_x = [1 size(entropy,2)];
tam_label = 6;
wi = 0.1;
hi = 0.07;
% tamaño de la grafica

labels = 0; % activo los numero de los ejes
fold = 1;
for s = 1:9
    for class = 1:2
        set(gcf,'position',[667   528   404   420])
        x1 = 0.02;     y1 = 0.06;     w = 0.95;      h = 0.92;
        subplot('position',[x1,y1,w,h]);
        load(['K:\ERD' filesep 'ERD_folds30_sub' num2str(s)]);
        ERDs = ERDsfilt_{fold}{class};
        limits_x = [0 1];
        limits_y = [0 1];
        main_topo_plots(ERDs,limits_x,limits_y,sel,tam_label,wi,hi,labels)
        saveas(gca,['G:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\figures' filesep 'ERDs_sub_' num2str(s) '_Class_' num2str(class)],'png')
        saveas(gca,['G:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\figures' filesep 'ERDs_sub_' num2str(s) '_Class_' num2str(class)],'epsc')
        close
    end
end
% save
% saveas(gcf,['head_s'],'epsc')