% Limpiar lo almacenado y la ventana de comandos.
clc; % close all
clearvars -except Disatfilt1_ Disatfilt2_
%% Dirección de la base de datos
% SUBJECTS_DIR = 'D:\BCI';
% clear
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
SUBJECTS_DIR5 = 'Particiones';                              % ubicación de las particiones
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

% Sujetos de BCICIV_2a
SS = 1:9;

% Sujetos de Giga
% SS = [37 32 12 18 42 34 3 7 35 33 21 2 4 39 29 43 28]; % UNO BUENO Y UNO MALO%%%% INDEXACDOS DE ACIERDO A CSP
% SS = [37 32 12 18 42 34 3 7];
% [37,15,7,1:6]; %6,14 [18:41]

% if strcmp(COHORT,'GIGASCIENCE_')   -   Mirar que sujetos no se utilizan.
%     SubInd = [50,14];
%     SS(SubInd) = [];
% end

%% paramaters definition
% Tiempo de referencia
tstart = 0.5;
tend = 1.5;
% definir parametros de filter bank
f_low  = 4;
f_high = 40;
Window = 4;
Ovrlap = 2;
filter_bank = [f_low:Ovrlap:f_high-Window;...
    f_low+Window:Ovrlap:f_high]';
orden_filter = 5; % orden del filtro.
labels = [1 2];   % labels utilizados.

% Numero de folds
Nfolds = 10;
% Carga de las particiones
load([SUBJECTS_DIR5 filesep 'cvNEW.mat'])

for s = SS
    reporte = [SUBJECTS_DIR4 filesep SUBJECTS{s} '.txt'];
    diary('on')
    diary(reporte)
    load([SUBJECTS_DIR filesep SUBJECTS{s} filesep 'eeg' filesep 'raw.mat'])
    y_ = y; a_ = ismember(y_,labels);
    y = y(:);
    ind = ismember(y,labels);
    y = y(ind);
    X = X(ind);
    X = cellfun(@(x) double(x) ,X,'UniformOutput',false);
    %     inicio = tic;
    %     % Filtro la señal.
    %     Xa = cell(size(filter_bank,1),1);
    %     diary(reporte)
    %     for b = 1:size(filter_bank,1)%Precompute all filters and trim
    %         Xa{b} = fcnfiltband(X, fs, filter_bank(b,:), 5);
    %         Xa{b} = cellfun(@(x) x(seg_start:seg_end,:),Xa{b},'UniformOutput',false);
    %     end
    % Creo los espacios.
    acc=nan(Nfolds,100,100);
    ks=nan(Nfolds,100,100);
    Xcp = cell(Nfolds,100);
    Xcp_ = cell(1,100);
    sfeats = cell(Nfolds,100);
    sfeats_ = cell(1,100);
    % Cargamos los ERDs de los sujetos.
    band = 'filter_ponderada';
    %     load(['I:\results rhy\ERDs_' num2str(s),'.mat'])
    load(['F:\New_ERD_abri_2019\ERD_folds_s',num2str(s)])
    % name =  'train'; % 'train'; %  'test'
    %   load(['G:\ERD_30\ERDs_Sub_',num2str(s),'_fold10filter_all',name,'.mat'],'r1')
    %     erd = r1; %
    pond = 0; % 0 es sin ponderar, 1 es ponderando
    aaaa =0;
    for len = 1.5%[2,1.5,1,0.5,0.2] % 0.5
        over = 0.1; % overlap
        fs = r1{1}{1}.fs;
        nsamples = numel(r1{1}{1}.t_plot);
        t = 1:round(len*fs*over):(nsamples-(len*fs));
        Npar1 = numel(t); %numero de ventanas
        %         clear r1
        %% Folds
        %         set(0,'DefaultFigureWindowStyle','docked')
        for Distan = 4%[1,3,4]
            %             figure
            for fold = 1:Nfolds
                clear Disa_
                % Distance computation- variable (Distan)
                % Distance: 1) inner product, 2) 1-inner product, 3) Euclidean,
                %4) Correlation, 5) 1-Correlation
                graf = 0; % graf = 1 si quieres graficar la mascara con la relevancia de los canales.
                distan_ = 1; %  distan_ = 1 si quiero calcular una distancia.
                nfreq = numel(erd{1}{1}.f_plot); 
                [Dis{s,fold}, erd_ord] = fncDistcomputation(erd{fold},len,over,nsamples,size(X{1},2),fs,[0,7],Distan,distan_,nfreq); % Dis (ch x freq x windows)
                % %                 freq = 1:17;
                %                 for ch = 1:size(Dis,1)
                %                     aa = 1;
                %                     for fre = fr1:fr2
                %                         for win = 1:size(Dis,3)
                %                             if fr1 == fr2
                %                                 Dis2(ch,win) = squeeze(Dis(ch,fre,win));
                %                             else
                %                                 Dis2(ch,aa,win) = Dis(ch,fre,win)./sum(Dis(ch,fr1:fr2,win));
                %                             end
                %                         end
                %                         aa = aa+1;
                %                     end
                %                 end
                
                %                 if fr1 == fr2
                %                     Dis_ = Dis2;
                %                 else
                %                     Dis_ = squeeze(sum(Dis2,2));
                %                 end
                
                % %por clase mirar en ritmos de eeg.
                % %                 if numel(freq) == 1
                % %                     for cl = 1:numel(labels)
                % %                         for ch = 1:size(erd_ord{cl},2)
                % %                             aa = 1;
                % %                             for fre = freq
                % %                                 for tim = 1:size(erd_ord{cl}{ch},2)
                % %                                     Disa_{cl}(ch,tim) = erd_ord{cl}{ch}(fre,tim);
                % %                                     aa = aa +1;
                % %                                 end
                % %                             end
                % %                         end
                % %                     end
                % %
                % %                     for cl = 1:numel(labels)
                % %                         for ch = 1:size(erd_ord{cl},2)
                % %                             tmp(ch,:) = Disa_{cl}(ch,:);
                % %                         end
                % %                         Disa_{cl} = tmp;
                % %                     end
                % %                     ERDsfilt_{fold} = Disa_;
                % %                     %                     t = 1:round(len*fs*over):(nsamples-(len*fs));
                % %                     %                     for cl = 1:numel(labels)
                % %                     %                         for ch = 1:size(Disa_{cl},1)
                % %                     %                             for tao = 1:length(t)
                % %                     %                                 Disaw_{cl}{ch,tao} = Disa_{cl}(ch,t(tao):(t(tao)+len*fs-1));
                % %                     %                             end
                % %                     %                         end
                % %                     %                     end
                % %                     %                     if aaaa == 1
                % %                     %                         Disa_ = Disatfilt1_;
                % %                     %                         for cl = 2%:numel(labels)
                % %                     %
                % %                     %                             figure;
                % %                     %                             for tao_ = [1,5,8,10,13,14,18,22,25]
                % %                     %                                 if tao_ == 1
                % %                     %                                     tao = tao_;
                % %                     %                                 else
                % %                     %                                     tao = t(tao_);
                % %                     %                                 end
                % %                     %
                % %                     %                                 load('G:\Dropbox\ERD\Codes\Topoplots\BCICIV_2a\electrodesBCICIV2a.mat')
                % %                     %                                 %load('BCICIV_2a\labels.mat')
                % %                     %                                 %load('BCICIV_2a\layout.mat')
                % %                     %                                 load('G:\Dropbox\ERD\Codes\Topoplots\HeadModel.mat')   % model of the head.
                % %                     %                                 Fil = 3; Col = 3;              % topoplot en subplots para cada uno de los sujetos.
                % %                     %                                 t_e = 70;                        % tamaño visual de los electrodos seleccionados.
                % %                     %                                 M1.xy = elec_pos(:,1:2);% posicion de los canales.
                % %                     %                                 M1.lab = Channels;        % nombre de los canales.
                % %                     %                                 % mask = reshape(rho,[17 22]) <= thresh; %
                % %                     %                                 sel = 1:22;%find(mask(band,:)==1); % Orden de los canales segun una relevancia.
                % %                     %                                 %  imagesc(1:1/250:7,1:22,Disa_{cl})
                % %                     %                                 rel = Disa_{cl}(:,tao);
                % %                     %                                 %                         rel = zeros(1,22);
                % %                     %                                 figure;
                % %                     %                                 fig = gca;
                % %                     %                                 %                 fig.Parent.PaperPosition = [2.91 4.15 2.69 2.69];
                % %                     %                                 %                             fig.Parent.OuterPosition = [1998 541 406 468];
                % %                     %                                 fig.Parent.OuterPosition = [2311 678 243 315];
                % %                     %                                 %                         rel = squeeze(er{cl}(:,t_(tim),fre));
                % %                     %                                 pos = M1.xy;
                % %                     %                                 label = M1.lab;
                % %                     %                                 % tam - tamaño de la cabeza y posición de los electrodos.
                % %                     %                                 tam = 1; tam2 = 0.96; tam3 = 0.98;
                % %                     %                                 for i=1:2
                % %                     %                                     pos(:,i) = tam2.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
                % %                     %                                 end
                % %                     %                                 xc = HeadModel(1,:);
                % %                     %                                 yc = HeadModel(2,:);
                % %                     %                                 %                         load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
                % %                     %                                 hold on
                % %                     %                                 % Topoplot
                % %                     %                                 x = pos(:,1);
                % %                     %                                 y = pos(:,2);
                % %                     %                                 tmp = [x,y,x*0]*rotz(2);
                % %                     %                                 x = tmp(:,1);
                % %                     %                                 y = tmp(:,2);
                % %                     %                                 pos(:,1) = x;
                % %                     %                                 pos(:,2) = y;
                % %                     %                                 GS = 500;
                % %                     %                                 xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                % %                     %                                 yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                % %                     %                                 [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                % %                     %                                 %% Creating data mask
                % %                     %                                 [TH R] = cart2pol(Xi,Yi);
                % %                     %                                 Zi(R>0.5) = NaN;
                % %                     %                                 deltax = xi(2)-xi(1); % length of grid entry
                % %                     %                                 deltay = yi(2)-yi(1); % length of grid entry
                % %                     %                                 %                             max(rel)
                % %                     %                                 h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                % %                     %                                 %                         shading interp % interpola colores.
                % %                     %                                 %                 scatter(tam3.*x(sel),(tam3+0.04).*y(sel),t_e,'b','filled')
                % %                     %                                 %                 text(pos(sel,1)*tam3-0.02,pos(sel,2)*tam3+0.02,label(sel),'Interpreter',...
                % %                     %                                 %                     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                % %                     %                                 %                     'FontSize',10);
                % %                     %                                 plot(tam*xc,tam*yc+0.001,'k','LineWidth',3)
                % %                     %                                 %                         colormap(erdcolormap)
                % %                     %                                 caxis([0 1])
                % %                     %                                 axis off
                % %                     %                                 hold off
                % %                     %                                 axis square
                % %                     %                                 axis image
                % %                     %                                 fs = 250;
                % %                     %                                 nsamples = 1751;
                % %                     %                                 over = 0.1;
                % %                     %                                 %                             colorbar
                % %                     %                                 t = 1:round(len*fs*over):(nsamples-(len*fs));
                % %                     %                                 title(['Sub ' num2str(s) ' cl ' num2str(cl) ' w ' num2str(t(tao_)/250)])
                % %                     %                                 axis xy
                % %                     %                                 %                             title(['Sub ' num2str(s) ' cl ' num2str(cl) ' w ' num2str(len)])
                % %                     %                             end
                % %                     %                         end
                % %                     %                     end
                % %                 else
                % %                     % por clase
                % %                     for cl = 1:numel(labels)
                % %                         if pond == 0
                % %                             %                             for ch = 1:size(erd_ord{cl},2)
                % %                             %                                 aa = 1;
                % %                             %                                 for fre = fr1:fr2
                % %                             %                                     rel_fre{cl}{ch}(aa)  = var(erd_ord{cl}{ch}(fre,:));
                % %                             %                                     aa = aa+1;
                % %                             %                                 end
                % %                             %                                 rel_fre{cl}{ch} = rel_fre{cl}{ch}./sum(rel_fre{cl}{ch});
                % %                             %                             end
                % %                             for ch= 1:size(erd_ord{cl},2)
                % %                                 aa = 1;
                % %                                 for fre = freq
                % %                                     Disa_{cl}{ch}(aa,:) = erd_ord{cl}{ch}(fre,:);
                % %                                     aa = aa+1;
                % %                                 end
                % %                             end
                % %                         else
                % %                             tmp = erd_ord{cl};
                % %                             tmpx = cellfun(@(x) x(:,0.5*fs:4.5*fs), tmp,'UniformOutput',false);
                % %                             tmpx = cell2mat(tmpx);
                % %                             pf = var(tmpx,[],2);% frequency relevence by channel
                % %                             pf = pf/max(pf);
                % %                             Dis_ = cell2mat(cellfun(@(x) sum((x),1), tmp,'UniformOutput',false)');
                % %                             for ch = 1:size(Dis_,1)
                % %                                 Disa_{cl}{ch}(1,:) = erd_ord{cl}{ch}(3,:);
                % %                                 Disa_{cl}{ch}(2,:) = Dis_(ch,:);
                % %                             end
                % %                             %                         figure; plot(tmpp')
                % %                         end
                % %                     end
                % %                     ERDsfilt_{fold} = Disa_;
                % % %                     for cl = 1:numel(labels)
                % % %                         for ch = 1:size(erd_ord{cl},2)
                % % %                             aa = 1;
                % % %                             for fre = fr1:fr2
                % % %                                 dats{cl}{ch}(aa,:) = erd_ord{cl}{ch}(fre,:).*rel_fre{cl}{ch}(aa); %datos escalados
                % % %                                 aa = aa +1;
                % % %                             end
                % % %                         end
                % % %                     end
                % % %                     for cl = 1:numel(labels)
                % % %                         for ch = 1:size(erd_ord{cl},2)
                % % %                             tmp(ch,:) = sum(dats{cl}{ch},1);
                % % %                         end
                % % %                         Disa_{cl} = abs(tmp)./max(max(abs(tmp)));
                % % %                     end
                % % %
                % % %                     t = 1:round(len*fs*over):(nsamples-(len*fs));
                % % %                     for cl = 1:numel(labels)
                % % %                         for ch = 1:size( Disa_{cl},1)
                % % %                             for tao = 1:length(t)
                % % %                                 Disaw_{cl}{ch,tao} = Disa_{cl}(ch,t(tao):(t(tao)+len*fs-1));
                % % %                             end
                % % %                         end
                % % %                     end
                % % %                     if aaaa == 1
                % % %                         Disa_=  Disatfilt2_;
                % % %                         for cl = 2%:numel(labels)
                % % %                             figure;
                % % %                             for tao_ = [1,5,8,10,13,14,18,22,25]
                % % %                                 if tao_ == 1
                % % %                                     tao = tao_;
                % % %                                 else
                % % %                                     tao = t(tao_);
                % % %                                 end
                % % %                                 load('G:\Dropbox\ERD\Codes\Topoplots\BCICIV_2a\electrodesBCICIV2a.mat')
                % % %                                 %load('BCICIV_2a\labels.mat')
                % % %                                 %load('BCICIV_2a\layout.mat')
                % % %                                 load('G:\Dropbox\ERD\Codes\Topoplots\HeadModel.mat')   % model of the head.
                % % %                                 Fil = 3; Col = 3;             % topoplot en subplots para cada uno de los sujetos.
                % % %                                 t_e = 70;               % tamaño visual de los electrodos seleccionados.
                % % %                                 M1.xy = elec_pos(:,1:2);% posicion de los canales.
                % % %                                 M1.lab = Channels;      % nombre de los canales.
                % % %                                 % mask = reshape(rho,[17 22]) <= thresh; %
                % % %                                 sel = 1:22;%find(mask(band,:)==1); % Orden de los canales segun una relevancia.
                % % %                                 rel = Disa_{cl}(:,tao);
                % % %                                 %                         rel = zeros(1,22);
                % % %                                 figure;
                % % %                                 fig = gca;
                % % %                                 %                 fig.Parent.PaperPosition = [2.91 4.15 2.69 2.69];
                % % %                                 %                             fig.Parent.OuterPosition = [1998 541 406 468];
                % % %                                 fig.Parent.OuterPosition = [2311 678 243 315];
                % % %                                 %                         rel = squeeze(er{cl}(:,t_(tim),fre));
                % % %                                 pos = M1.xy;
                % % %                                 label = M1.lab;
                % % %                                 % tam - tamaño de la cabeza y posición de los electrodos.
                % % %                                 tam = 1; tam2 = 0.96; tam3 = 0.98;
                % % %                                 for i=1:2
                % % %                                     pos(:,i) = tam2.*((pos(:,i)-min(pos(:,i)))/(range(pos(:,i)))-0.5);
                % % %                                 end
                % % %                                 xc = HeadModel(1,:);
                % % %                                 yc = HeadModel(2,:);
                % % %                                 %                         load('G:\Dropbox\ERD\results_ERDfc_subjects\Mapas de colores\erdscolormap.mat')
                % % %                                 hold on
                % % %                                 % Topoplot
                % % %                                 x = pos(:,1);
                % % %                                 y = pos(:,2);
                % % %                                 tmp = [x,y,x*0]*rotz(2);
                % % %                                 x = tmp(:,1);
                % % %                                 y = tmp(:,2);
                % % %                                 pos(:,1) = x;
                % % %                                 pos(:,2) = y;
                % % %                                 GS = 500;
                % % %                                 xi = linspace(min(x)-0.05, max(x)+0.05, GS);       % x-axis for interpolation (row vector)
                % % %                                 yi = linspace(min(y)-0.05, max(y)+0.05, GS);       % y-axis for interpolation (row vector)
                % % %                                 [Xi, Yi, Zi] = griddata([x' xc], [y' yc], [rel(:);zeros(numel(xc),1)], xi', yi,'v4'); % interpolate the topographic data
                % % %                                 %% Creating data mask
                % % %                                 [TH R] = cart2pol(Xi,Yi);
                % % %                                 Zi(R>0.5) = NaN;
                % % %                                 deltax = xi(2)-xi(1); % length of grid entry
                % % %                                 deltay = yi(2)-yi(1); % length of grid entry
                % % %                                 %                             max(rel)
                % % %                                 h = surf(Xi-deltax/2, Yi-deltay/2+0.004, zeros(size(Zi)), Zi,'EdgeColor', 'none', 'FaceColor', 'flat');hold on
                % % %                                 %                         shading interp % interpola colores.
                % % %                                 %                 scatter(tam3.*x(sel),(tam3+0.04).*y(sel),t_e,'b','filled')
                % % %                                 %                 text(pos(sel,1)*tam3-0.02,pos(sel,2)*tam3+0.02,label(sel),'Interpreter',...
                % % %                                 %                     'latex','ButtonDownFcn',{@lineCallback,database},'Color','black',...
                % % %                                 %                     'FontSize',10);
                % % %                                 plot(tam*xc,tam*yc+0.001,'k','LineWidth',3)
                % % %                                 %                         colormap(erdcolormap)
                % % %                                 %                             caxis([0 1])
                % % %                                 axis off
                % % %                                 hold off
                % % %                                 axis square
                % % %                                 axis image
                % % %                                 fs = 250;
                % % %                                 nsamples = 1751;
                % % %                                 over = 0.1;
                % % %                                 t = 1:round(len*fs*over):(nsamples-(len*fs));
                % % %                                 title(['Sub ' num2str(s) ' cl ' num2str(cl) ' w ' num2str(t(tao_)/250)])
                % % %                                 axis xy
                % % %                                 %                             title(['Sub ' num2str(s) ' cl ' num2str(cl) ' w ' num2str(len)])
                % % %                             end
                % % %                         end
                % % %                     end
                [~,pos1]=min(abs(t./250-2));
                [~,pos2]=min(abs(t./250-4.5));
                rho = mean(Dis{s,fold}(:,:,pos1:pos2),3)'; %
                rho = rho(:);
                threshold = zeros(1,100);
                threshold = linspace(min(rho),max(rho),100);
                [Pi_f{s,fold},mask_i{s,fold}] = fncmask(rho,threshold,filter_bank,s,X,graf,fold);
                % %
                % %                 end
                % %
                % %                 %Calculo de la mascara según las ventanas que se quieren
                % %                 %analizar.
                % %                 %                 mass = 0;
                % %                 %                 if mass == 1
                % %                 %             rho = mean(Dis{s,fold}(:,:,11:23),3)'; %
                % %                 %             rho = rho(:);
                % %                 %             threshold = zeros(1,100);
                % %                 %             threshold = linspace(min(rho),max(rho),100);
                % %                 %             [Pi_f{s,fold},mask_i{s,fold}] = fncmask(rho,threshold,filter_bank,s,X,graf,fold);
                % %                 %             %                     % tr_ind = c{s,fold}.training; tr_ind = tr_ind(a_);
                % %                 %                     % ts_ind = c{s,fold}.test; ts_ind = ts_ind(a_);
                % %                 %                     % banda = cell(Nfolds,100);
                % %                 %                     % for u = 1:numel(threshold)
                % %                 %                     %
                % %                 %                     % end
                % %                 %                 elseif mass == 2
                % %                 %                     % pendiente
                % %                 %                 end
            end
            %
            %             if numel(freq) == 1
            % %                 save(['S:\Mi unidad\ERDs_filtrados',filesep,'ERDs_Sub_',num2str(s),'_fold',num2str(fold),band],'ERDsfilt_')
            %                 save(['D:\Google Drive\ERDs_filtrados',filesep,'ERDs_Sub_',num2str(s),'_fold',num2str(fold),band],'ERDsfilt_')
            %             else
            % %                 save(['S:\Mi unidad\ERDs_filtrados',filesep,'ERDs_Sub_',num2str(s),'_fold',num2str(fold),band],'ERDsfilt_')
            %                 save(['D:\Google Drive\ERDs_filtrados',filesep,'ERDs_Sub_',num2str(s),'_fold',num2str(fold),band],'ERDsfilt_')
            %             end
            %
            %             graf = 0;
            %             if Distan == 1; Met = 'innerproduct'; elseif  Distan == 2; Met = '1-innerproduct';
            %             elseif  Distan == 3; Met = 'Euclidean'; elseif  Distan == 4; Met = 'Correlation';
            %             else Met = '1-Correlation '; end
            %             if graf == 1
            %                 suptitle(['Subject: ',num2str(s),' Method:',Met]);  colorbar('Position',[0.93 0.1 0.01 0.8]);
            %                 saveas(gcf,['Graficas',filesep,'Subj_',num2str(s),'_',Met,'_',num2str(len*1000),'beta'],'png')
            %                 close
        end
        
    end
end
% end


% % clear Di_cor
% % for cl = 1
% %     aa = 1;
% %     for tao_ = 1:1751%[1,5,8,10,13,14,18,22,25]
% %         %         if tao_ == 1
% %         tao = tao_;
% %         %         else
% %         %             tao = t(tao_);
% %         %         end
% %         Di_cor{cl}(aa) = abs(corr(sum(Disatfilt2_{cl}(:,tao),2),sum(Disatfilt1_{cl}(:,tao),2)));
% %         aa = aa+1;
% %
% %     end
% % end
% % hold on; plot(0:1/250:7,Di_cor{cl})
% % hold on; plot(0:1/250:7,erd_ord{2}{8}(2,:),'v')
% % title(['Sujeto 8 clase ',num2str(cl)])
% % legend('Rw,nw','Rw,f','Rnw,f','ERD')
% Disatfilt1_, Disatfilt2_

