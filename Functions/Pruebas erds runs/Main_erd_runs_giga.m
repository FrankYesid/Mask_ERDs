%% Template for using the EEG analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment: Estimation ERDs in runs.
% ERD: evente-related desynchronization.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019 Signal Processing and Recognition Group.
% Universidad Nacional de Colombia.
% L.F. Velasquez-Martinez
% F.Y. Zapata-Castaño.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: run BIOSIG for functions of load dataset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc %
warning off
%% path toolbox biosig
run('F:\Dropbox\ERD\results_ERDfc_subjects\pruebas_mask_ERDs\Toolbox\Biosig_ERD\biosig_installer.m')
addpath(genpath('F:\Dropbox\ERD\results_ERDfc_subjects\pruebas_mask_ERDs\Toolbox\Biosig_ERD\biosig\t310_ERDSMaps'))
%% Load data
% name of database
name = 'GigaScience';
% Path database
SUBJECTS_DIR = 'F:\Gigascience\data';
% format of data.
data = '.mat'; % download in Gigascience

%% Parameters
%% ------------- ERDs -----------------
lambda = 0;
alpha  = 0.01;     % alpha of significance in ERDs.
t      = [0, 0, 7];% Start point, time resolution and end point within a trial (in s)
% <1x3>. If the second value is 0, the time resolution corresponds to 1/fs.
ref    = [0.5,1.5];% reference of estimation ERDs.
cue    = 2;        % cue of dataset.

%% Select GigaScience.
SS     = flip([1:28,30:33,35:52]);  % Numer of Subjects in
Class_ = 1:2; % Select classes of database.

COHORT = 's';
SUBJECTS = dir([SUBJECTS_DIR filesep '*' COHORT '*']);
SUBJECTS = struct2cell(SUBJECTS);
SUBJECTS = SUBJECTS(1,:)';
% elimina los sujetos de Gigascience
if strcmp(COHORT,'s') && numel(SS)==numel(1:52)
    SubInd = [29,34]; % eliminamos los sujetos con poca cantidad de ensayos necesarios para el calculo del ERDs.
    SS(SubInd) = [];
end

%% all runs
load Result\trial_sub_giga.mat
% one run consists in 40 trials (20 for each of the two possible classes).
channels_select = [13,21,31,48,50,58];
Runs   = cell(1,numel(SS));
indx_  = cell(numel(SS),numel(Class_));
num_trials = cell(numel(SS),numel(Class_));
% 2 runs por cada clase.
% for ss = SS
%     for class = 1:2
%         temp = in(ss,class);
%         tem  = floor(temp/2);
%         dife = temp-tem; 
%         if dife-tem == 0
%             num_trials{ss,class}{1} = temp/2;% numero de intentos por cada clase prueba.
%             num_trials{ss,class}{2} = temp/2;% numero de intentos por cada clase prueba.
%         else
%             num_trials{ss,class}{1} = dife;% numero de intentos por cada clase prueba.
%             num_trials{ss,class}{2} = tem;% numero de intentos por cada clase prueba.
%         end
%     end
% end
% con 30 trials por cada run.
Ntrials = 30;
for ss = SS
    for class = 1:2
        temp = in(ss,class);
        if mod(temp,Ntrials) == 0
            rn = temp/Ntrials;
            tem = cell(rn,1);
            tem_= cellfun(@(x) Ntrials,tem,'UniformOutput',false);
        else
            rn = floor(temp/Ntrials)+1;
            tem = cell(rn,1);
            tem_= cellfun(@(x) Ntrials,tem,'UniformOutput',false);
            tem_{numel(tem_)} = mod(temp,Ntrials);
        end 
        num_trials{ss,class} = tem_;
    end
end
% calculo de los ERDs por runs
for ss = SS
    load([SUBJECTS_DIR filesep SUBJECTS{ss}])
    trials = cell(numel(Class_),2);
    r1 = cell(1,numel(Class_));
    erd = cell(1,numel(Class_));
    for clas = Class_              % Clases Utilizadas.
        num_tr       = num_trials{ss,clas};
        method       = 'bp';       % 'bp' or 'fft'
        refmethod    = 'classic';  % 'classic' or 'trial'
        bad_trails   = eeg.bad_trial_indices.bad_trial_idx_mi{clas};
        h.Classlabel = ones(numel(find(eeg.imagery_event==1)),1)*clas;
        h.SampleRate = eeg.srate;
        h.TRIG       = find(eeg.imagery_event==1)';
        indd         = h.Classlabel; indd(bad_trails) = 0;
        h.Classlabel = h.Classlabel(logical(indd));
        inde = indd(logical(indd));
        indx_        = zeros(numel(find(eeg.imagery_event==1)),1);
        ve_run = [];
        for ve = 1:numel(num_trials{ss,clas})
            ve_run = [ve_run;ones(num_tr{ve},1)*ve];%;ones(num_tr,1).*3;...
        end
        %             ones(num_tr,1).*4;ones(num_tr,1).*5];
        for runs = 1:numel(num_trials{ss,clas})%2%in(ss,clas)/num_tr
            trials{clas,runs} = num_trials{ss,clas}{runs}; %
            tic
            if clas == 1
                S = double(eeg.imagery_left'); % clase de la imaginación de la mano izquierda.
                S(isnan(S)) = 0;  % Coloca en cero los datos NAN que se encuentren en le señal.
                % Montaje de la base de datos.
                montage = '64ch'; % montaje de los canales que se van a utilizar.
                channels = 1:64;  % número de canales que se van a utilizar.
                % Output parameters:
                % lap         ... Laplacian filter matrix.
                % plot_index  ... Indices for plotting the montage.
                % n_rows      ... Number of rows of the montage.
                % n_cols      ... Number of columns of the montage.
                [lap, plot_index, n_rows, n_cols] = getMontage(montage);
                mont = plot_index;% Como cuadrar la información relacionada.
                % Multiplicamos la señal por la matriz de filtrado laplaciano.
                s   = S(:,channels)*lap;
                s   = s(:,channels_select);
%                 cv  = ones(1,numel(find(eeg.imagery_event==1)));
%                 cv  = cv(1:numel(find(eeg.imagery_event==1))-numel(bad_trails));
                t_t = numel(find(eeg.imagery_event==1))-numel(bad_trails);
%                 a   = 1;
%                 for i = 1:numel(h.Classlabel)
%                     if inde(i) == 1
%                         indx_(i) = cv(a);
%                         a = a+1;
%                     else
%                         indx_(i) = 1;
%                     end;
%                 end;
                indx = inde & ismember(ve_run,runs);
            elseif clas == 2
                S = double(eeg.imagery_right');  % clase de la imaginación de la mano derecha.
                S(isnan(S)) = 0;    % Coloca en cero los datos NAN que se encuentren en le señal.
                % Montaje de la base de datos.
                montage = '64ch'; % montaje de los canales que se van a utilizar.
                channels = 1:64;   % número de canales que se van a utilizar.
                % Output parameters:
                %   lap        ... Laplacian filter matrix.
                %   plot_index ... Indices for plotting the montage.
                %   n_rows     ... Number of rows of the montage.
                %   n_cols     ... Number of columns of the montage.
                [lap, plot_index, n_rows, n_cols] = getMontage(montage);
                mont = plot_index;  % Como cuadrar la información relacionada.
                s   = S(:,channels)*lap;
                s   = s(:,channels_select);
                % Multiplicamos la señal por la matriz de filtrado laplaciano.
                cv = ones(1,numel(find(eeg.imagery_event==1))-numel(bad_trails));
%                 a = 1;
%                 for i = 1:numel(h.Classlabel)
%                     if inde(i) == 2
%                         indx_(i) = cv(a);
%                         a = a+1;
%                     else
%                         indx_(i) = 0;
%                     end;
%                 end;
                indx = inde & ismember(ve_run,runs);
            end
            %% ERDS maps
            r1{clas}{runs} = calcErdsMap(indx,s, h, [-2, 0, 5], [19, 19], 'method', method, 'class', clas, 'ref', [1.5,2], 'f_bandwidths', 22, 'f_steps', 22, 'sig', 'boxcox', 'lambda', 0, 'alpha', 0.01, 'heading', name, 'montage', mont, 'cue', 2, 'refmethod', refmethod);
            r1{clas}{runs}.heading = ['Subject: ' num2str(ss)];
            erd{clas}{runs} = plotErdsMap2(r1{clas}{runs},ss,clas,10);
            fprintf(['Subject...',num2str(ss),' Class...',num2str(clas),' Run...',num2str(runs),' Time...',num2str(toc),'\n'])
        end
    end
    ERDs_runs_con = cell(1,2);
    for clas = Class_
        for run = 1:numel(erd{clas})
            for ch = 1:numel(erd{clas}{run}.ERDS)
                ERDs_runs_con{clas}{run}{ch} =  erd{clas}{run}.ERDS{ch}.erds';
            end
        end
    end
    ERDs_runs_sin = cell(1,2);
    for clas = Class_
        for run = 1:numel(erd{clas})
            for ch = 1:numel(erd{clas}{run}.ERDS)
                ERDs_runs_sin{clas}{run}{ch} =  r1{clas}{run}.ERDS{ch}.erds';
            end
        end
    end
    save(['F:\Sub_' num2str(ss) '_ERD_runs_',num2str(Ntrials),'_trials_2.mat'],'ERDs_runs_con','ERDs_runs_sin','trials')
end