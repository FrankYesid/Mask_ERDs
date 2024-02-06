% SFBCSP + MIL + for over several window sizes for more instances
clc; clear; close all;
warning off;
%--------------------------------------------------------------------------
% Database
subjects_dir = 'H:\OpenNeuro';
% subjects_dir = '/home/dfcollazosh/Matlab_Tests/Pruebas_MIL_LASSO_OpenNeuro/openNeuro/';
%--------------------------------------------------------------------------
% Functions
addpath(genpath([ subjects_dir '/functions/']));
filesave = 'G:/OpenNeuro/resultadosERP/';

%--------------------------------------------------------------------------
% cohort
cohort   = 'sub-';
subjects = dir([subjects_dir filesep '*' cohort '*']);
subjects = struct2cell(subjects);
subjects = subjects(1,:)';
name = {'aud'; 'vis'};
%--------------------------------------------------------------------------

ch{1} = {8 24};
ch{2} = {5 30};
ch{3} = {29 30 34 36};
ch{4} = {23 24 30};
ch{5} = {8 30};
ch{6} = {9 24 30};
ch{7} = {30};
ch{8} = {28 29 30 33 36};
ch{9} = {29 30 33};
ch{10} = {24 29 33};
ch{11} = {16 24 30};
ch{12} = {30};
ch{13} = {30};
ch{14} = {30};
ch{15} = {8 16 28 34 37};
ch{16} = {16 24 30 28 37};
ch{17} = {7 33};

% subjects
ss = (1:length(subjects));
% chan_locs

for s = ss
    %------------------------------------------------------------------
    % search grid and labels to use
    labels = [0 1];
    %------------------------------------------------------------------
    clearvars -except ch X_all tw1 conteo1 Total_Acc conteo2 tw_min tw_max ...
        twi tws tw s TVW TVWa mat_aport ss mmacc mmaccTr Names Acc table ...
        labels Twindow lags param filter_bank experiment_name subjects ...
        cohort subjects_dir filesave name
    display(subjects{s})
    load('D:\Dropbox\ERD\Codes\Topoplots\Channels 34 channels\M_loc.mat');
    
    tic
    %--------------------------------------------------------------
    % load data and labels
    X = []; Y = [];
    for task = 1
        for run = 1:3
            load([subjects_dir filesep subjects{s} filesep 'EEG' filesep ['task00' num2str(task) '_run00' num2str(run)] filesep 'EEG_noGA.mat']);
            fs               = 1000;
            XT               = data_noGA(1:43,:);
            %% Filter notch, power line
            Xfil             = fcnfiltband_matrix(XT, fs, 1, 5, 3); % high-pass
            Xfil             = fcnfiltband_matrix(Xfil, fs, 60, 35, 2); % notch
            Xfil             = fcnfiltband_matrix(Xfil, fs, 120, 5, 2); % notch
            Xfil             = fcnfiltband_matrix(Xfil, fs, 100, 5, 1); % low-pass
            XT = Xfil;
            
            %% Filter BCG
            Xfil         = fcnfiltband_matrix(Xfil, fs, 4, 5, 1); % low-pass
            [W, y]        = pca(Xfil');
            tmp = zeros(43,1);
            tmp([1 2]) = 1;
            IC = diag(tmp);
            X_ = W\IC*y';
            XT1              = XT - X_;
            
            %% Re-Reference
            X_den            = shortestpath(XT1, cell2mat(ch{s}));
            
            %% Build trials
            number_of_trials = 125;
            trial_time       = load([subjects_dir filesep subjects{s} ...
                filesep 'EEG' filesep ['task00' num2str(task) '_run00' num2str(run)] filesep 'trial_time.txt']);
            t                = linspace(0,size(X_den,2)/fs,size(X_den,2));
            tao              = 1000;
            XT               = trials_construction(X_den,trial_time,number_of_trials,t,tao)';
            X                = [X;XT];
            YT               = load([subjects_dir filesep subjects{s} filesep 'EEG' filesep ['task00' num2str(task) '_run00' num2str(run)] filesep 'labels.txt']);
            Y                = [Y;YT];
        end
        
        tmp = X(Y == 1);
        tar = permute(reshape(cell2mat(tmp),  34, size(tmp,1), 1000), [2 1 3]);
        tmp = X(Y == 0);
        non_tar = permute(reshape(cell2mat(tmp),  34, size(tmp,1), 1000), [2 1 3]);
        
        %                 savefile = [filesave num2str(s) '/epochs_' name{task} '/tar.mat'];
        %                 save(savefile, 'tar', '-v7.3');
        %
        %                 savefile = [filesave num2str(s) '/epochs_' name{task} '/non_tar.mat'];
        %                 save(savefile, 'non_tar', '-v7.3');
        %                 eeglab; close
        if task == 1
           save(['D:\Subject',num2str(s)],'tar','non_tar')
        else
           save(['D:\Subject',num2str(s),'_2'],'tar','non_tar')
        end
        channels = [1:34];
        EEG = struct_EEGLAB(tar,channels);
        typeplot = 1;
        figure;
        % Fz = 7;
        erpimage( mean(EEG.data([25], :),1), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), 'PZ', 10, 1 ,'yerplabel','\muV','erp','on','cbar','on','topo', { [25] EEG.chanlocs EEG.chaninfo } );
        figure; pop_erpimage(EEG,typeplot)%,typeplot,[],trial_time)
        figure;pop_timtopo(EEG, [-200 800]);
        
        X = []; Y = []; tmp = [];
        saveas(gca,['task_',num2str(task)],'epsc')
    end
    clearvars data_noGA t
    
end
