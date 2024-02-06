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
%% Load data
% name of database
name = 'BCICIV_2a';
% Path database
path_name = 'K:\Database\BCI Competition\BCICIV_2a\Training\A0';
% format of data. 
data = '.gdf'; % download in BCI competition IV.
% cvpartition of the trials.
load partition\cvNEW.mat
%% Parameters
Nfolds = 1;    % Numer of folds

%% ------------- Filter bank -------------
f_low  = 4;    % frecuency minimal.
f_high = 40;   % frecuency maxim.
f_bandwidths = 4;    % bandwith.
Ovrlap = 2;    % overlape.
filter_bank = [f_low:Ovrlap:f_high-f_bandwidths;...
    f_low+f_bandwidths:Ovrlap:f_high]';

%% ------------- ERDs -----------------
lambda = 0;
alpha  = 0.01;     % alpha of significance in ERDs.
t      = [0, 0, 7];% Start point, time resolution and end point within a trial (in s) 
                   % <1x3>. If the second value is 0, the time resolution corresponds to 1/fs.
ref    = [0.5,1.5];% reference of estimation ERDs.
cue    = 2;        % cue of dataset.
%%  BCI Competition IV dataset 2a.
SS     = 1:9; % Numer of Subjects in
Class_ = 1:2; % Select classes of database.
% Method   ... 1 Select trials of train.
%          ... 2 Select trials of test.
%          ... 3 Select All trals.
Method_ = 3;

%% all runs 
% one run consists in 48 trials (12 for each of the four possible classes).
runs = 1:6;
ve_run = [ones(48,1);ones(48,1).*2;ones(48,1).*3;ones(48,1).*4;ones(48,1).*5;ones(48,1).*6];
Runs = cell(1,numel(SS));
indx_ = cell(numel(SS),numel(Class_));
for ss = SS
    if Method_ < 3        
        erd_sin = cell(1,Nfolds);
        erd_con =cell(1,Nfolds);
        for fold = 1:Nfolds         % folds.
            for clas = Class_       % Class.
                classes = [1 2 3 4];% [] ... All classes in database.
                method = 'bp';      % 'bp' or 'fft'.
                refmethod = 'classic';  % 'classic' or 'trial'.
                channel = 0;       % 0 ... All channels
                % function of biosig - load subjects of dataset.
                [S, h] = sload([path_name num2str(ss) 'T' data], channel, 'OVERFLOWDETECTION:OFF');
                S(isnan(S)) = 0;
                montage = '22ch';  % montage electrodes of database.
                channels = 1:22;   % numer of channels.
                %% Output parameters:
                %   lap            ... Laplacian filter matrix.
                %   plot_index     ... Indices for plotting the montage.
                %   n_rows         ... Number of rows of the montage.
                %   n_cols         ... Number of columns of the montage.
                [lap, plot_index, n_rows, n_cols] = getMontage(montage);
                mont = plot_index;
                % Laplacian filter.
                s = S(:,channels)*lap;
                % Delete of artifacts.
                indd = h.ArtifactSelection==0 & ismember(h.Classlabel, classes);
                indx = ones(1,numel(h.TRIG(indd)))';
                h.Classlabel = h.Classlabel(indd);
                h.TRIG = h.TRIG(indd);
                % select of trails respect of cvpartition.
                if Method_ == 1     % Train.
                    indx = indx.*c{ss,fold}.training.*ismember(h.Classlabel,clas);
                elseif Method_ == 2 % Test.
                    indx = indx.*c{ss,fold}.test.*ismember(h.Classlabel,clas);
                end                
                fs = h.SampleRate;
                cue = 1;
                erd_sin{fold}{clas} = calcErdsMap(indx,s, h, t, [filter_bank(1,1)+2, filter_bank(end,1)+2],...
                    'method', method, 'class', clas, 'ref', ref,'f_bandwidths', f_bandwidths, ...
                    'f_steps', Ovrlap, 'sig', 'boxcox', 'lambda', lambda,'alpha', alpha, 'heading', name,...
                    'montage', mont, 'cue', cue,'refmethod', refmethod, 'submean', true);
                erd_sin{fold}{clas}.heading = ['Subject: ' num2str(ss)]; % Name of subject.
                % Quantification ERD/ERS with significance.
                erd_con{fold}{clas} = plotErdsMap2(erd_sin{fold}{clas},ss,clas,channels);
                toc
            end
        end
    else
        erd_sin = cell(1,numel(Class_));   
        erd_con =cell(1,numel(Class_));   
        indx = [];
        for clas = Class_ 
            tic
            classes = [1 2 3 4]; % [] ... All classes - Clases para los indices.
            method = 'bp';     % 'bp' or 'fft'.
            refmethod = 'classic';  % 'classic' or 'trial'.
            % Carga base de datos.
            channel = 0;      % 0 ... All channels
            [S, h] = sload([path_name num2str(ss) 'T' data], channel, 'OVERFLOWDETECTION:OFF');
            S(isnan(S)) = 0;    % Coloca en cero los datos NAN que se encuentren en le señal.
            % Montaje de la base de datos.
            montage = '22ch'; % montaje de los canales que se van a utilizar.
            channels = 1:22;   % número de canales que se van a utilizar.
            % Output parameters:
            %   lap            ... Laplacian filter matrix.
            %   plot_index ... Indices for plotting the montage.
            %   n_rows     ... Number of rows of the montage.
            %   n_cols      ... Number of columns of the montage.
            [lap, plot_index, n_rows, n_cols] = getMontage(montage);
            mont = plot_index;  % Como cuadrar la información relacionada.
            % Multiplicamos la señal por la matriz de filtrado laplaciano.
            s = S(:,channels)*lap;
            % Organizamos la base de datos de los trials que se quieren utilizar de la base de datos.
            indd = h.ArtifactSelection==0 & ismember(h.Classlabel, classes);
            h.Classlabel = h.Classlabel(indd);
            h.TRIG = h.TRIG(indd);
            fs = h.SampleRate;
%             for r = runs
%                 indx = ismember(h.Classlabel,[1 2]) & ismember(ve_run(indd),r);
%                 indx_{s,clas} = sum(indx);
%                 % Quantification ERD/ERS.
%                 erd_sin{clas}{r} = calcErdsMap(indx,s, h, t, [filter_bank(1,1)+2, filter_bank(end,1)+2],...
%                     'method', method, 'class', clas, 'ref', ref,'f_bandwidths', f_bandwidths, ...
%                     'f_steps', Ovrlap, 'sig', 'boxcox', 'lambda', lambda,'alpha', alpha, 'heading', name,...
%                     'montage', mont, 'cue', cue,'refmethod', refmethod, 'submean', true);
%                 erd_sin{clas}.heading = ['Subject: ' num2str(ss)];
%                 % Quantification ERD/ERS with significance.
%                 erd_con{clas}{r} = plotErdsMap2(erd_sin{clas}{r},ss,clas,channels);
%                 fprintf('Subject ')
%             end
        end
    end
%     if Method_ == 1
%         name1 = 'train';
%     elseif Method_ == 2
%         name1 = 'test';
%     else
%         name1 = 'all';
%     end
%     fprintf(['Quantification Event-related Desynchronization Subject: ',num2str(ss),'Type: ',name1,' in ',name,'\n'])
%     ERDs_con = cell (1,Nfolds);
%     ERDs_sin = cell (1,Nfolds);
%     if numel(erd_con) == Nfolds
%         for fold = 1:Nfolds
%             for cl = Class_
%                 for ch = channels
%                     ERDs_con{fold}{cl}{ch} = erd_con{fold}{cl}.ERDS{ch}.erds;
%                 end
%             end
%         end
%         for fold = 1:Nfolds
%             for cl = Class_
%                 for ch = channels
%                     ERDs_sin{fold}{cl}{ch} = erd_sin{fold}{cl}.ERDS{ch}.erds;
%                 end
%             end
%         end
%         save(['Result\ERD_folds',num2str(Nfold),'_sub' num2str(ss) nam '.mat'],'ERDs_con','ERDs_sin')
%     else
%         for fold = 1:Nfolds
%             for cl = Class_
%                 for ch = channels
%                     ERDs_con{fold}{cl}{ch} = erd_con{fold}{cl}.ERDS{ch}.erds;
%                 end
%             end
%         end
%         for fold = 1:Nfolds
%             for cl = Class_
%                 for ch = channels
%                     ERDs_sin{fold}{cl}{ch} = erd_sin{fold}{cl}.ERDS{ch}.erds;
%                 end
%             end
%         end
%         save(['Result\ERD_sub' num2str(ss) nam 'class.mat'],'ERDs_con','ERDs_sin')
%     end
    Runs{ss} = ve_run(indd);
    clear erd_sin erd_con indd
end
save('Result\Runs_trials.mat','Runs')
