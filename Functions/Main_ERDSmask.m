% Limpiar lo almacenado y la ventana de comandos.
clear; close all; clc

%% Dirección de la base de datos
% SUBJECTS_DIR = 'D:\BCI';

% ---------------------------- LF ----------------------------
% SUBJECTS_DIR = 'D:\BCI';                                  % Base de datos
% SUBJECTS_DIR2 = 'D:\Luisa\Dropbox\ERD\results_ERDfc_subjects\BCI'; Almacena resultados de acc en dropbox.
% SUBJECTS_DIR3 = 'D:\Luisa\Dropbox\ERD\ERD_2019\BCI';      % Almacena graficas de ERD en Dropbox.
% SUBJECTS_DIR4 = 'D:\Luisa\Dropbox\ERD\Codigos Corriendo'; % Almacena las los reportes.
% % SUBJECTS_DIR4 = 'D:\BCI';                                 % particiones se encuentran en la misma carpte de BCI.
% SUBJECTS_DIR5 = 'Particiones';                            % particiones - se encuentra en Dropbox.
% SUBJECTS_DIR6 = 'Graficas';                                 % Guardar Graficas de acc u otras.

% ---------------------------- fy1 ---------------------------
SUBJECTS_DIR = 'F:\BCI';                                    % Carga de Base de datos.
SUBJECTS_DIR2 = 'G:\Dropbox\ERD\results_ERDfc_subjects\BCI';% Almacena resultados de acc en dropbox.
SUBJECTS_DIR3 = 'G:\Dropbox\ERD\ERD_2019\BCI';              % Almacena graficas de ERD en Dropbox.
SUBJECTS_DIR4 = 'G:\Dropbox\ERD\Codigos Corriendo';         % Almacena las los reportes.
SUBJECTS_DIR5 = 'Particiones\BCICIV_2a';                              % ubicación de las particiones
SUBJECTS_DIR6 = 'Graficas';                                 % Guardar Graficas de acc u otras.

%% Direccion del fold de las funciones
% ---------------------------- LF ----------------------------
% addpath(genpath('D:\Luisa\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('D:\Luisa\Dropbox\ERD\results_ERDfc_subjects'))
% ---------------------------- LF-pc -------------------------
% addpath(genpath('C:\Users\lfvelasquezm\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('C:\Users\lfvelasquezm\Desktop\frank\functions'))

% ---------------------------- fy1 ---------------------------
addpath(genpath('G:\Dropbox\ERD\Codes\TP\Matlab_wang\csp\CSP_fun\functions'))
% addpath(genpath('G:\Dropbox\ERD\results_ERDfc_subjects'))
% ---------------------------- fy2 ---------------------------
% addpath(genpath('C:\Users\frany\Dropbox\Event-related\Codes\TP\Matlab_wang\csp\CSP_fun\functions'));
% addpath(genpath('C:\Users\frany\Dropbox\Event-related\results_ERDfc_subjects'))

%% DataBase
% % BCIIII_4a_
% % BCICIV_2a_
%  GIGASCIENCE_

COHORT = 'BCICIV_2a_';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';

% %% grilla de busqueda para lambda
param = linspace(0,0.9,100);
experiment_name = mfilename;

if strcmp(COHORT,'BCICIV_2a_')
    % Sujetos de BCICIV_2a
    SS = 1:9;
else
    % Sujetos de Giga
    %    SS = 1:52;
    SS = [37 32 12 18 42 34 3 7 35 33 21 2 4 39 29 43 28]; % UNO BUENO Y UNO MALO%%%% INDEXACDOS DE ACIERDO A CSP
    %    SS = [37 32 12 18 42 34 3 7];
    %    [37,15,7,1:6]; %6,14 [18:41]
    %    if strcmp(COHORT,'GIGASCIENCE_')   % -   Mirar que sujetos no se utilizan y los quita de la lista.
    %        SubInd = [50,14];
    %        SS(SubInd) = [];
    %    end
end

%% paramaters definition
% Tiempo de referencia - para el calculo del ERDs.
tstart = 0.5;
tend = 1.5;
% definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
orden_filter = 5;  % orden del filtro.
labels = [1 2];     % labels utilizados.
% Numero de folds
Nfolds = 30;
% Carga de las particiones
% load([SUBJECTS_DIR5 filesep 'cvNEW.mat'])
load(['G:\Dropbox\ERD\results_ERDfc_subjects\cv_new_giga.mat'])

Distancia_ = cell(1,numel(SS));
for s = SS
    reporte = [SUBJECTS_DIR4 filesep SUBJECTS{s} '.txt'];
    % Guarda el reporte del codigo que se ejecuta
    diary('on')
    diary(reporte)
    % Carga de los datos del sujeto.
    load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    y_ = y; a_ = ismember(y_,labels);
    y = y(:);
    ind = ismember(y,labels);
    y = y(ind);
    X = X(ind);
    X = cellfun(@(x) double(x) ,X,'UniformOutput',false);
    inicio = tic;
    % Filter bank.
    Xa = cell(size(filter_bank,1),1);
    diary(reporte)    
    for b = 1:size(filter_bank,1)%Precompute all filters and trim
        Xa{b} = fcnfiltband(X, fs, filter_bank(b,:), 5);
        Xa{b} = cellfun(@(x) x(seg_start:seg_end,:),Xa{b},'UniformOutput',false);
    end % filter bank    
    % Creo los espacios de las variables.
    acc=nan(Nfolds,100,100);
    ks=nan(Nfolds,100,100);
    Xcp = cell(Nfolds,100);
    Xcp_ = cell(1,100);
    sfeats = cell(Nfolds,100);
    sfeats_ = cell(1,100);
    threshold = zeros(1,100);
    % Cargamos los ERDs del sujeto, con significancia y con aplicación de filtro Laplaciano.
    load(['F:\New_ERD_abri_2019\ERD_folds' num2str(s) '_2.mat'])
    for len = 2 % 0.5, 1, 1.5, 2.
        over = 0.1;
        fs = r1{1}{1}.fs;
        nsamples = numel(r1{1}{1}.t_plot);
        t = 1:round(len*fs*over):(nsamples-(len*fs));
        Npar1 = numel(t); %numero de ventanas
        clear r1
        %% Folds
        %         set(0,'DefaultFigureWindowStyle','docked')
        for Distan = 1:5 % Distance computation
            % Distance: 1) inner product, 2) 1-inner product, 3) Euclidean,
            %4) Correlation, 5) 1-Correlation
            mask_i = cell(1,Nfolds);
            Dis = cell(1,Nfolds);
            for fold = 1:Nfolds                
                graf = 0;       % graf = 1 si quieres graficar la mascara con la relevancia de los canales.
                distan_ = 1; %  distan_ = 1 si quiero calcular una distancia.
                MI = 1;         % 1 - si se mira solo la ventana de MI, 0 - si se realiza varias ventanas con un translape definido anteriormente.
                nfreq =  numel(erd{1}{1}.f_plot);
                [Dis{fold}, ~] = fncDistcomputation(erd{fold}, len,over,nsamples,size(X{1},2),fs,[0,7],Distan,distan_,nfreq,1); % Dis (ch x freq x windows)
                %Calculo de la mascara según las ventanas que se quieren
                %analizar.
                if numel(size(Dis{fold})) < 3
                    rho = mean(Dis{fold},3)'; % importante saber si se calcula el promedio en las ventanas de MI o en todas las ventanas de tiempo.
                else
                    rho = Dis{fold};
                end
                rho = rho(:);
                threshold = linspace(min(rho),max(rho),100);
                [Pi_f,mask_i{fold}] = fncmask(rho,threshold,filter_bank,s,X,graf,fold);
                % tr_ind = c{s,fold}.training; tr_ind = tr_ind(a_);
                % ts_ind = c{s,fold}.test; ts_ind = ts_ind(a_);
                % banda = cell(Nfolds,100);
                % for u = 1:numel(threshold)
                %                
                %         end
                %     end
                % end
                %         rho = abs(erd_mean(:))./max(erd_mean(:));
                %         threshold = zeros(1,100);
                %         threshold = linspace(min(rho),max(rho),100);
                %         % threshold(end)=1;
                %         if sum(isnan(rho))>1; rho(isnan(rho)==1)=1;  end %para los nan detectados
                %         % tr_ind   = cv.training(fold); tr_ind = tr_ind(ind);
                %         % ts_ind   = cv.test(fold); ts_ind = ts_ind(ind);
                %         tr_ind = c{s,fold}.training; tr_ind = tr_ind(a_);
                %         ts_ind = c{s,fold}.test; ts_ind = ts_ind(a_);
                %         banda = cell(Nfolds,100);
                %         for u = 1:numel(threshold)
                %             tic
                %             mask = reshape(rho,[size(filter_bank,1) size(X{1},2)]) <= threshold(u);
                %             rho,[size(filter_bank,1) size(X{1},2)]
                %             % mask = reshape(p,[33,22])<= threshold(u);
                %             SelBand = sum(mask,2);
                %             valQ = floor(SelBand/2);
                %             selected = SelBand>=6;%6 seleccuionde canales en esa frecuencia
                %             band = find(selected==1);
                % %             banda{fold,u}=band;
                %             if sum(selected) == 0
                %                 continue
                %             end
                %             Xc = cell(1,numel(band));
                %             for b = 1:numel(band)
                %                 [~,chan]=find(mask(band(b),:)>0);
                %                 if numel(chan) < 6 %1
                %                     continue
                %                 end
                %                 C = cell2mat(reshape(cellfun(@(x)(cov(x(:,chan))/trace(cov(x(:,chan)))),Xa{band(b)},'UniformOutput',false),[1 1 numel(Xa{b})]));
                %                 W = csp_feats(C(:,:,tr_ind),y(tr_ind),'train','Q',3);%floor(numel(chan)/2)
                %                 Xc{b} = csp_feats(C,W,'test');
                %             end
                %             Xc = cell2mat(Xc);
                %             Xcp{fold,u} = Xc;
                %             clear C W
                %             % Lasso
                %             target = mapminmax(y(tr_ind)')';
                %             B = lasso(Xc(tr_ind,:),target,'Lambda',param);
                %             selected_feats = abs(B) > eps;
                %             sfeats{fold,u} = selected_feats;
                %             %
                %             parfor l=1:numel(param)
                %                 Xcc = Xc(:,selected_feats(:,l));
                %                 if size(Xcc,2)<2
                %                     continue
                %                 end
                %                 mdl = fitcdiscr(Xcc(tr_ind,:),y(tr_ind)); %LDA
                %                 acc(fold,u,l) = mean(mdl.predict(Xcc(ts_ind,:))==reshape(y(ts_ind),[sum(ts_ind) 1]));
                %                 % Confusion Matrix
                %                 tar_pred = mdl.predict(Xcc(ts_ind,:)); %tar_pred(tar_pred==1)=0; tar_pred(tar_pred==2)=1;
                %                 tar_true = reshape(y(ts_ind),[sum(ts_ind) 1]); %tar_true(tar_true==1)=0; tar_true(tar_true==2)=1;
                %                 conM = confusionmat(tar_true,tar_pred);
                %                 ks(fold,u,l) = kappa(conM);
                %                 % plotconfusion(tar_true',tar_pred');
                %             end % lambda
                %             fprintf(['Threshold...' num2str(u) '...' num2str(toc) '\n'])
                %             %             [fold,u]
                %         end % signficance
                %         fprintf(['Sujeto: ' SUBJECTS{s} ' Fold...' num2str(fold) '\n'])
                %         % toc
                %         diary(reporte)
            end % folds
            Distancia_{s}{Distan} = Dis;
        end % distancia
        %
        %     act = squeeze(mean(acc,1));
        %     [dato,indp] = max(act(:));
        %     actstd = squeeze(std(acc,1)); actstd = actstd(:); actstd = actstd(indp);
        %     [u_opt,l_opt]=ind2sub(size(act),indp);
        %     table(1,:) = [threshold(u_opt),param(l_opt),dato*100,actstd*100];
        %     save([SUBJECTS_DIR2 ...
        %          filesep SUBJECTS{s} filesep experiment_name 'Results_vLF.mat'],'acc','table','ks');
        %     fprintf([' ...acc: ' num2str(dato*100,'%02.1f') ' std: ' num2str(actstd*100,'%02.1f')...
        %         ' ...time: ' num2str(toc-inicio) '\n']);
        %     diary(reporte)
        %     %limpia los datos que se vuelven a calcular.
        %     clear acc table Xcp sfeats ks
        %
        %     % Seccion de la validacion de los sujetos
        %     % Cargamos los datos de evaluación de cada sujeto.
        %     load(['F:\BCI\BCI2a_Evaluation' filesep 'BCI_s0' num2str(s) 'E.mat'])
        %     X_ = cell(size(BCI_data.BCI_data,3),1);
        %     for ch = 1:size(BCI_data.BCI_data,1)
        %         for tr = 1:size(BCI_data.BCI_data,3)
        %             X_{tr}(:,ch) = BCI_data.BCI_data(ch,:,tr);
        %         end
        %     end
        %     %Filtrado de la señal de evaluacion.
        %     Xa_ = cell(size(filter_bank,1),1);
        %     for b = 1:size(filter_bank,1)%Precompute all filters and trim
        %         Xa_{b} = fcnfiltband(X_, fs, filter_bank(b,:), 5);
        %         Xa_{b} = cellfun(@(x) x(seg_start:seg_end,:),Xa_{b},'UniformOutput',false);
        %     end
        %     y_ = labels.label; ind_ = ismember(y_,[1,2]); y_ = y_(ind_);
        %     tr_ind = ones(numel(c{s,1}.training),1); tr_ind = tr_ind(a_);
        %     ts_ind = ones(numel(y_),1); %ts_ind = ts_ind(a_);
        %
        %     u = u_opt; % threshold optimo
        %     inicio2 = tic;
        %     % como organizar un solo rho? -----------------------------------------
        %     mask = reshape(rho,[size(filter_bank,1) size(X{1},2)]) <= threshold(u);
        %     % mask = reshape(p,[33,22])<= threshold(u);
        %     SelBand = sum(mask,2); valQ = floor(SelBand/2);
        %     selected = SelBand>=6; %6 seleccuionde canales en esa frecuencia
        %     band = find(selected==1); banda{u}=band;
        %     Xc = cell(1,numel(band)); Xc_ = cell(1,numel(band));
        %     for b = 1:numel(band)
        %         [~,chan]=find(mask(band(b),:)>0);
        %         if numel(chan) < 6 % selecciono 6 canales 1
        %             continue
        %         end
        %         % Para train
        %         C = cell2mat(reshape(cellfun(@(x)(cov(x(:,chan))/trace(cov(x(:,chan)))),Xa{band(b)},'UniformOutput',false),[1 1 numel(Xa{b})]));
        %         W = csp_feats(C(:,:,tr_ind),y(tr_ind),'train','Q',3);%floor(numel(chan)/2)
        %         Xc{b} = csp_feats(C,W,'test');
        %         % Para test
        %         C_ = cell2mat(reshape(cellfun(@(x)(cov(x(:,chan))/trace(cov(x(:,chan)))),Xa_{band(b)},'UniformOutput',false),[1 1 numel(Xa_{b})]));
        %         W_ = csp_feats(C_(:,:,ts_ind),y_(ts_ind),'train','Q',3);%floor(numel(chan)/2)
        %         Xc_{b} = csp_feats(C_,W_,'test');
        %     end
        %     Xc = cell2mat(Xc); Xcp{u} = Xc; clear C W
        %     Xc_ = cell2mat(Xc_); Xcp_{u} = Xc_; clear C_ W_
        %     % Lasso
        %     target = mapminmax(y(tr_ind)')';
        %     target_ = mapminmax(y_(ts_ind)')';
        %     B = lasso(Xc(tr_ind,:),target,'Lambda',param);
        %     B_ = lasso(Xc_(ts_ind,:),target_,'Lambda',param);
        %     selected_feats = abs(B) > eps;
        %     selected_feats_ = abs(B_) > eps;
        %
        %     sfeats{u} = selected_feats;
        %     sfeats_{u} = selected_feats_;
        %     %Clasificación
        %     l=l_opt; % lambda optimo
        %     Xcc = Xc(:,selected_feats(:,l));
        %     Xcc_ = Xc_(:,selected_feats_(:,l));
        %     mdl = fitcdiscr(Xcc(tr_ind,:),y(tr_ind));      %LDA
        %     acc_ = mean(mdl.predict(Xcc_(ts_ind,:))==reshape(y_(ts_ind),[sum(ts_ind) 1]));
        %     % Confusion Matrix
        %     tar_pred_ = mdl.predict(Xcc_(ts_ind,:));         %tar_pred(tar_pred==1)=0; tar_pred(tar_pred==2)=1;
        %     tar_true_ = reshape(y_(ts_ind),[sum(ts_ind) 1]); %tar_true(tar_true==1)=0; tar_true(tar_true==2)=1;
        %     conM_ = confusionmat(tar_true_,tar_pred_);
        %     ks_ = kappa(conM_);
        %     % plotconfusion(tar_true',tar_pred');
        %     % end % lambda
        %     %     fprintf(['Threshold...' num2str(u) '...' num2str(toc) '\n'])
        %     %   [fold,u]
        %     % end % signficance
        %     table(1,:) = [threshold(u_opt),param(l_opt),dato*100,actstd*100];
        %     save([SUBJECTS_DIR2 ...
        %          filesep SUBJECTS{s} filesep experiment_name 'Results_vLF.mat'],'acc','table','ks');
        %     fprintf([' ...acc: ' num2str(dato*100,'%02.1f') ' std: ' num2str(actstd*100,'%02.1f')...
        %         ' ...time: ' num2str(toc-inicio) '\n']);
        %     diary(reporte)
        %     diary(reporte)
        %     diary('off')
        %     clearvars -except filename filter_bank SS Nfolds
    end % tamaño de la ventana.    
end % Sujeto 

for s = [2,8]
    for Distan = [1,3,4]
        figure;
        a = 1;
        for  fold = 1:10%Nfolds
            subplot(2,5,a)
            imagesc(Distancia_{s}{Distan}{fold})
            a = a+1;
        end
        if Distan == 1
            suptitle('Prod. Punto')
            name = 'Prod. Punto';
        elseif Distan == 3
            suptitle('Euclidea')
            name = 'Euclidea';
        elseif Distan == 4
            suptitle('Corr')
            name = 'Corr';
        end
        saveas(gcf,['G:\Dropbox\ERD\Main_ERDS\Graficas\Mask\Dis',name,'_sub',num2str(s)],'epsc')
        saveas(gcf,['G:\Dropbox\ERD\Main_ERDS\Graficas\Mask\Dis',name,'_sub',num2str(s)],'png')
    end
end

% % Graficar los aciertos de cada sujeto. tanto para train como para test
% % Hacer la tabla de los aciertos - en la tabla ya existente de acc.
% %
% %                 save(['C:\Users\Lufe\Desktop\Prueba_dif',filesep,'Sub_',num2str(s),'_w',num2str(len*1000),'_Met',num2str(Distan),'_fold',num2str(fold)],'Dis_')
% %             end
% %             graf = 0;
% %             if Distan == 1; Met = 'innerproduct'; elseif  Distan == 2; Met = '1-innerproduct';
% %             elseif  Distan == 3; Met = 'Euclidean'; elseif  Distan == 4; Met = 'Correlation';
% %             else Met = '1-Correlation '; end
% %             if graf == 1
% %                 suptitle(['Subject: ',num2str(s),' Method:',Met]);  colorbar('Position',[0.93 0.1 0.01 0.8]);
% %                 saveas(gcf,['Graficas',filesep,'Subj_',num2str(s),'_',Met,'_',num2str(len*1000),'alpha'],'png')
% %                 close
% %             end
