%% Graficar los ERD/ERS
% Tenemos en cuenta que se grafican los datos en EEG.
%% BCI_competition_IV dataset 2a
%% Gráfica de los 3 canales centrales c3, cz, c4
SS = [8,2]; load('E:\copia disco D\Escritorio\DRop\pruebas_mask_ERDs\Mapas de colores\erdscolormap.mat')

nam = '100';
for s = SS
    %     figure(s)
    load(['I:\New_ERD_abri_2019\Train_folds\Estructura\ERD_folds' num2str(s) '_2.mat'])
    savea = 1;  % 1... para guardar la image n, 0 ... no guarda la imagen.
    Dir_save = 'G:\Dropbox\ERD\ERD_2019\Informe_general\Figs\ERDs_maps\BCICIV_2a\channels_center\'; % ubicación donde queremos guardar la gráfica de los sujetos.
    %     if s == 3 || s == 1 || s == 7 || s == 6 || s == 2 || s == 4  || s == 9   % selecciona sujetos que no tienen label y
    labels = 1;
    %     else
    %         labels = 1;
    %     end
    %     if s == 8 || s == 9 || s == 3 || s == 1 || s == 5 || s == 7 % selecciona sujetos que no tienen
    labels2 = 1;
    %     else
    %         labels2 = 0;
    %     end
    load('maxmin.mat')
    r = erd;
    for fold = 1
        for cl = 1:2
            a = 1;
            for ch = [8,10,12]
                %                 r{fold}{cl}.ERDS{ch}.erds = (-1+ ((erd{fold}{cl}.ERDS{ch}.erds-min(Mini{s}(fold,cl,a))).*(1.5-(-1))/(max(Maxi{s}(fold,cl,a))-min(Mini{s}(fold,cl,a)))))*100;
                %                 for fr = 1:17;r{fold}{cl}.ERDS{ch}.erds(:,fr) = r{fold}{cl}.ERDS{ch}.erds(:,fr)-mean(r{fold}{cl}.ERDS{ch}.erds(:,fr)) ; end;

                r{fold}{cl}.ERDS{ch}.erds = ((erd{fold}{cl}.ERDS{ch}.erds)./abs(min(Mini{s}(:))));
                a = a+1;
            end
        end
    end
    Graficas_ERDs_channels(r,erdcolormap,labels,labels2,s,savea,Dir_save);
    %     saveas(gca,['C:\Users\frany\Desktop\Documento_tesis\Figures\Figs_cap1\ERD_maps\Subjects_',num2str(s)],'png')
    matlab2tikz(['D:\Documento_tesis\Figures\Figs_cap1\ERD_maps\Subjects2_',num2str(s),'.tikz'])
    pause(0.01)
    %     close
    clear r1 erd
end

%% BCI_competition_IV dataset 2a_ Class 3 and 4
%% Gráfica de los 3 canales centrales c3, cz, c4
SS = [8,1,2]; load('H:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
close all
for s = SS
    %     figure(s)
    %     load(['F:\New_ERD_abri_2019\Train_folds\ERD_folds' num2str(s) '.mat'])
    %     load(['F:\New_ERD_abri_2019\ERD_folds30_sub',num2str(s),'_class3_4.mat'])
    load(['L:\ERD\ERD_folds30_sub',num2str(s),'.mat'])
    savea = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'H:\Dropbox\ERD\Main_ERDS\Graficas\'; % ubicación donde queremos guardar la gráfica de los sujetos.
    if s == 1 || s == 2    % selecciona sujetos que no tienen label y
        labels = 0;
    else
        labels = 1;
    end
    if s == 8 || s == 1 || s == 2 % selecciona sujetos que no tienen
        labels2 = 1;
    else
        labels2 = 0;
    end
    %     load('maxmin.mat')
    %     r = erd;
    %     for fold = 1
    %         for cl = 1:2
    %             a = 1;
    %             for ch = [8,12]
    % %                 r{fold}{cl}.ERDS{ch}.erds = (-1+ ((erd{fold}{cl}.ERDS{ch}.erds-min(Mini{s}(fold,cl,a))).*(1.5-(-1))/(max(Maxi{s}(fold,cl,a))-min(Mini{s}(fold,cl,a)))))*100;
    % %                 for fr = 1:17;r{fold}{cl}.ERDS{ch}.erds(:,fr) = r{fold}{cl}.ERDS{ch}.erds(:,fr)-mean(r{fold}{cl}.ERDS{ch}.erds(:,fr)) ; end;
    %
    %                 r{fold}{cl}.ERDS{ch}.erds = ((erd{fold}{cl}.ERDS{ch}.erds)./abs(min(Mini{s}(:))));
    %                 a = a+1;
    %             end
    %         end
    %     end
    Graficas_ERDs_channels2(ERDsfilt_,erdcolormap,labels,labels2,s,savea,Dir_save);
    pause(0.01)
    close
    clear ERDsfilt_ ERDsfilt__o
end

%% Gráfica del sujeto en los 22 canales.
SS = [2,8];
% load('D:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
nam = '100';
for s = SS
    load(['I:\New_ERD_abri_2019\Train_folds\Estructura\ERD_folds' num2str(s) '_2.mat'])
    save = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    for clas = 1:2
        Dir_save = 'G:\Dropbox\ERD\Main_ERDS\Graficas\Laplacian\'; % ubicación donde queremos guardar la gráfica de los sujetos.
        Graficas_ERDs_channels_all(erd,erdcolormap,s,save,Dir_save,clas)
        close
    end
    clear r1 erd
end

%% Gráfica del sujeto en los 22 canales_version2.
SS = 1:9; load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    load(['F:\New_ERD_abri_2019\ERD_folds30_sub',num2str(s),'_class3_4.mat'])
    %     load(['F:\New_ERD_abri_2019\Train_folds\ERD_folds' num2str(s) '_2.mat'])
    save = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    for clas = 1:2
        Dir_save = 'G:\Dropbox\ERD\ERD_2019\Informe_general\Figs\ERDs_maps\BCICIV_2a\all_channels\';
        %         Dir_save = 'G:\Dropbox\ERD\Main_ERDS\Graficas\Laplacian\'; % ubicación donde queremos guardar la gráfica de los sujetos.
        Graficas_ERDs_channels_all2(ERDsfilt_,erdcolormap,s,save,Dir_save,clas)
        close
    end
    clear ERDsfilt_ ERDsfilt__o
end

%% Gráfica de los sujetos en la variación de los trials
SS = 8; load('D:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    load(['F:\New_ERD_abri_2019\Train_folds\ERD_folds' num2str(s) '_2.mat'])
    save = 0;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'G:\Dropbox\ERD\Main_ERDS\Graficas\Laplacian\'; % ubicación donde queremos guardar la gráfica de los sujetos.
    Graficas_ERDs_channels_trials(erd,erdcolormap,s,save,Dir_save)
    close
    clear r1 erd
end

%% Gráfica de las mascaras
SS = 1:9; load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    load(['F:\New_ERD_abri_2019\Train_folds\ERD_folds' num2str(s) '_2.mat'])
    save = 0;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'G:\Dropbox\ERD\Main_ERDS\Graficas\Laplacian\'; % ubicación donde queremos guardar la gráfica de los sujetos.
    Grafica_mask(erd,erdcolormap,s,save,Dir_save)
    close
    clear r1 erd
end

%% GigaScience
%% Gráfica de 64 canales de GigaScience
%% definir parametros de filter bank
SS = 14;%[1:28,30:33,35:52];
% SS = [17];
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
% load('D:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
Dir_save = 'D:\Dropbox\ERD\ERD_2019\Informe\Figs\ERDs_maps_all\Giga'; % ubicación donde queremos guardar la gráfica de los sujetos.
for s = SS
    load(['I:\ERD_giga\MI\ERD_folds1_s' num2str(s) '.mat'])
    save = 0;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    for clas = 1:2
        Graficas_ERDs_channels_all_giga(ERDsfilt_,erdcolormap,s,save,Dir_save,clas)
        close
    end
    clear r1 erd
end

%% Gráfica de los 3 canales centrales c3, cz, c4
% SS = [1:28,30:33,35:52];
SS = [46];% [1,27,38,43];
% SS = 11;
load('H:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    %     figure(s)
    load(['L:\ERD_giga\MI\30_folds_training\ERD_folds_s' num2str(s) '.mat'])
    savea = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'H:\Dropbox\ERD\Main_ERDS\Graficas\ERDs_giga\ERDs_maps_3ch\Giga'; % ubicación donde queremos guardar la gráfica de los sujetos.
    if s == 14 || s == 46 || s == 38 || s == 27 % selecciona sujetos que no tienen label y
        labels = 0;
    else
        labels = 1;
    end
    if s == 43 || s == 14 || s == 46 % selecciona sujetos que tienen eje x
        labels2 = 1;
    else
        labels2 = 0;
    end
    %     load('maxmin.mat')
    %     r = erd;
    for fold = 1
        for cl = 1:2
            a = 1;
            for ch = [13,14,48,50]
                %                 r{fold}{cl}.ERDS{ch}.erds = (-1+ ((erd{fold}{cl}.ERDS{ch}.erds-min(Mini{s}(fold,cl,a))).*(1.5-(-1))/(max(Maxi{s}(fold,cl,a))-min(Mini{s}(fold,cl,a)))))*100;
                %                 for fr = 1:17;r{fold}{cl}.ERDS{ch}.erds(:,fr) = r{fold}{cl}.ERDS{ch}.erds(:,fr)-mean(r{fold}{cl}.ERDS{ch}.erds(:,fr)) ; end;
                r{fold}{cl}.ERDS{ch}.erds = (ERDsfilt_{fold}{cl}{ch});%./abs(min(Mini{s}(:))));
                %                 r{fold}{cl}{ch} = (ERDsfilt_{fold}{cl}{ch});%./abs(min(Mini{s}(:))));
                a = a+1;
            end
        end
    end
    Graficas_ERDs_channels_giga2(r,erdcolormap,labels,labels2,s,savea,Dir_save);
    pause(0.01)
    close
    clear r ERDsfilt_
end

%% Gráfica de los canales centrales
SS = [1:28,30:33,35:52];
% SS = 11;
load('H:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    %     figure(s)
    load(['L:\ERD_giga\MI\30_folds_training\ERD_folds_s' num2str(s) '.mat'])
    savea = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'H:\Dropbox\ERD\Main_ERDS\Graficas\ERDs_giga\ERDs_maps_zone_MI\Giga'; % ubicación donde queremos guardar la gráfica de los sujetos.
    %     if s == 3 || s == 1 || s == 7 || s == 5 || s == 6 || s == 2 % selecciona sujetos que no tienen label y
    %         labels = 0;
    %     else
    labels = 1;
    %     end
    %     if s == 8 || s == 9 || s == 3 || s == 1 || s == 5 || s == 7 % selecciona sujetos que tienen eje x
    labels2 = 1;
    %     else
    %         labels2 = 0;
    %     end
    %     load('maxmin.mat')
    %     r = erd;
    for fold = 1
        for cl = 1:2
            a = 1;
            for ch = [14,13,12,48,49,50,51]
                %                 r{fold}{cl}.ERDS{ch}.erds = (-1+ ((erd{fold}{cl}.ERDS{ch}.erds-min(Mini{s}(fold,cl,a))).*(1.5-(-1))/(max(Maxi{s}(fold,cl,a))-min(Mini{s}(fold,cl,a)))))*100;
                %                 for fr = 1:17;r{fold}{cl}.ERDS{ch}.erds(:,fr) = r{fold}{cl}.ERDS{ch}.erds(:,fr)-mean(r{fold}{cl}.ERDS{ch}.erds(:,fr)) ; end;
                r{fold}{cl}.ERDS{ch}.erds = (ERDsfilt_{fold}{cl}{ch});%./abs(min(Mini{s}(:))));
                %                 r{fold}{cl}{ch} = (ERDsfilt_{fold}{cl}{ch});%./abs(min(Mini{s}(:))));
                a = a+1;
            end
        end
    end
    Graficas_ERDs_channels_giga3(r,erdcolormap,labels,labels2,s,savea,Dir_save);
    %         pause(0.01)
    close
    clear r ERDsfilt_
end


%% Gráfica de los 2 canales centrales c3, c4 all giga
SS = [1:28,30:33,35:52];
load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
nam = '100';
for s = SS
    %     figure(s)
    %     load(['K:\ERD_giga\Santiago\ERDs\ERD_30folds_s' num2str(s) '.mat'])
    load(['L:\ERD_giga\Santiago\ERDs',filesep,'ERD_30folds_s',num2str(s),'.mat'])
    savea = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
    Dir_save = 'H:\Dropbox\ERD\ERD_2019\Informe_general\Figs\ERDs_maps\Giga\channels_c3_c4\'; % ubicación donde queremos guardar la gráfica de los sujetos.
    %     if s == 3 || s == 1 || s == 7 || s == 5 || s == 6 || s == 2 % selecciona sujetos que no tienen label y
    %         labels = 0;
    %     else
    labels = 1;
    %     end
    %     if s == 8 || s == 9 || s == 3 || s == 1 || s == 5 || s == 7 % selecciona sujetos que tienen eje x
    labels2 = 1;
    %     else
    %         labels2 = 0;
    %     end
    %     load('maxmin.mat')
    %     r = erd;
    for fold = 1
        for cl = 1:2
            a = 1;
            for ch = [1,2]
                %                 r{fold}{cl}.ERDS{ch}.erds = (-1+ ((erd{fold}{cl}.ERDS{ch}.erds-min(Mini{s}(fold,cl,a))).*(1.5-(-1))/(max(Maxi{s}(fold,cl,a))-min(Mini{s}(fold,cl,a)))))*100;
                %                 for fr = 1:17;r{fold}{cl}.ERDS{ch}.erds(:,fr) = r{fold}{cl}.ERDS{ch}.erds(:,fr)-mean(r{fold}{cl}.ERDS{ch}.erds(:,fr)) ; end;
                r{fold}{cl}.ERDS{a}.erds = (ERDsfilt_{fold}{cl}{ch});%./abs(min(Mini{s}(:))));
                a = a+1;
            end
        end
    end
    Graficas_ERDs_channels_giga(r,ERDsfilt_,erdcolormap,labels,labels2,s,savea,Dir_save);
    pause(0.01)
    close
    clear r1 erd
end

%% Other
%% grafica del ERD modificando frecuencias
load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
nam = '100';
figure;fig = imagesc(0:1/250:7,1:20,ERDsfilt_{1}{1}{8}, [-1, 1.5]);
colormap(erdcolormap)
fig.Parent.YTick = 1:20;
fig.Parent.YTickLabel = num2str([r1{1, 1}{1, 1}.f_low;r1{1, 1}{1, 1}.f_up]');
axis xy; title(['ERD sub ',num2str(ss),' class ',num2str(1),' C3'])
figure;fig = imagesc(0:1/250:7,1:20,ERDsfilt_{1}{2}{8}, [-1, 1.5]);
colormap(erdcolormap)
fig.Parent.YTick = 1:20;
fig.Parent.YTickLabel = num2str([r1{1, 1}{1, 1}.f_low;r1{1, 1}{1, 1}.f_up]');
axis xy; title(['ERD sub ',num2str(ss),' class ',num2str(2),' C3'])
figure;fig = imagesc(0:1/250:7,1:20,ERDsfilt_{1}{1}{12}, [-1, 1.5]);
colormap(erdcolormap)
fig.Parent.YTick = 1:20;
fig.Parent.YTickLabel = num2str([r1{1, 1}{1, 1}.f_low;r1{1, 1}{1, 1}.f_up]');
axis xy; title(['ERD sub ',num2str(ss),' class ',num2str(1),' C4'])
figure;fig = imagesc(0:1/250:7,1:20,ERDsfilt_{1}{2}{12}, [-1, 1.5]);
colormap(erdcolormap)
fig.Parent.YTick = 1:20;
fig.Parent.YTickLabel = num2str([r1{1, 1}{1, 1}.f_low;r1{1, 1}{1, 1}.f_up]');
axis xy; title(['ERD sub ',num2str(ss),' class ',num2str(2),' C4'])

%% Gráfica de conectividad canales centrales
% for cl = 1:2
%     load(['E:\Cx_s_all\Cx_s8_c' num2str(cl) '.mat'])
%     for freq = 1:size(Cx_{cl},3)
%         for time = 1:size(Cx_{cl},4)
%                 CX{cl}(:,freq,time) = strengths_und(threshold_proportional(abs(X_1(:,:,freq,time)),0.3));
%         end
%     end
% end
close all
thr = 0.1;
SS = [2,8];
for s = SS
    load(['I:\cx_fun_2a\Cx_wpli_BCI_2a_all_time_0_7_sub' num2str(s) '_folds_1.mat'])
    for cl = 1:2
        for freq = 1:size(Cx_{cl},3)
            for time = 1:size(Cx_{cl},4)
                CX{cl}(:,freq,time) = strengths_und(threshold_proportional(abs(Cx_{cl}(:,:,freq,time)),thr));
            end
        end


        load('D:\Dropbox\ERD\results_ERDfc_subjects\ERDs_en topoplots\erdscolormap.mat')
        savea = 0;  % 1... para guardar la image n, 0 ... no guarda la imagen.
        Dir_save = 'E:\ERD_2019\Informe_general\Figs\Cx_maps\BCICIV_2a\channels_center\'; % ubicación donde queremos guardar la gráfica de los sujetos.
        %     if s == 3 || s == 1 || s == 7 || s == 6 || s == 2 || s == 4  || s == 9   % selecciona sujetos que no tienen label y
        labels = 1;
        %     else
        %         labels = 1;
        %     end
        %     if s == 8 || s == 9 || s == 3 || s == 1 || s == 5 || s == 7 % selecciona sujetos que no tienen
        labels2 = 1;
        %     else
        %         labels2 = 0;
        %     end
        % load('maxmin.mat')
        fold = 1;

        %     for cl = 1:2
        for ch = 1:22
            r{fold}{cl}.ERDS{ch}.erds = squeeze(CX{cl}(ch,:,:))';
        end
        r{fold}{cl}.classes = cl;
        r{fold}{cl}.class = cl;
        r{fold}{cl}.refmethod = 'all';
        r{fold}{cl}.t_plot = 0:1/250:7;
        % ------------- Filter bank -------------
        f_low  = 4;    % frecuencia mínima.
        f_high = 40;  % frecuencia máxima.
        Window = 4; % Tamaño de la banda de frecuencia.
        Ovrlap = 2;    % Traslape de las bandas.
        filter_bank = [f_low:Ovrlap:f_high-Window;...
            f_low+Window:Ovrlap:f_high]';
        r{fold}{cl}.f_plot = mean(filter_bank,2);
        r{fold}{cl}.f_low = filter_bank(:,1)';
        r{fold}{cl}.f_up = filter_bank(:,2)';
        %     end
        colorm = flip(colormap("hot"));
        %     [maxi,mini]=Graficas_Cx_channels(r,colorm,labels,labels2,s,savea,Dir_save);
        [maxi,mini]= Graficas_Cx_channels_all(r,colorm,s,savea,Dir_save,cl);
        close all

        load(['I:\cx_fun_2a\Cx_wpli_BCI_2a_all_time_0_7_sub' num2str(s) '_folds_1.mat'])
        %     for cl = 1:2
        for freq = 1:size(Cx_{cl},3)
            for time = 1:size(Cx_{cl},4)
                CX{cl}(:,freq,time) = strengths_und(threshold_proportional(abs(Cx_{cl}(:,:,freq,time)),thr));
            end
        end
        %     end
        savea = 1;  % 1... para guardar la imagen, 0 ... no guarda la imagen.
        fold = 1;
        %     for cl = 1:2
        for ch = 1:22
            r{fold}{cl}.ERDS{ch}.erds = squeeze(CX{cl}(ch,:,:))'/max(maxi(:));
        end
        r{fold}{cl}.classes = cl;
        r{fold}{cl}.class = cl;
        r{fold}{cl}.refmethod = 'all';
        r{fold}{cl}.t_plot = 0:1/250:7;
        % ------------- Filter bank -------------
        f_low  = 4;    % frecuencia mínima.
        f_high = 40;  % frecuencia máxima.
        Window = 4; % Tamaño de la banda de frecuencia.
        Ovrlap = 2;    % Traslape de las bandas.
        filter_bank = [f_low:Ovrlap:f_high-Window;...
            f_low+Window:Ovrlap:f_high]';
        r{fold}{cl}.f_plot = mean(filter_bank,2);
        r{fold}{cl}.f_low = filter_bank(:,1)';
        r{fold}{cl}.f_up = filter_bank(:,2)';
        %     end
        colorm = flip(colormap("hot"));
        %     [maxi,mini]=Graficas_Cx_channels(r,colorm,labels,labels2,s,savea,Dir_save);
        [maxi,mini]= Graficas_Cx_channels_all(r,colorm,s,savea,Dir_save,cl);
        close
    end
end