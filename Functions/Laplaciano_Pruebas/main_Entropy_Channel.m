%% Pruebas Correcion Articulo
close all; clear; clc
addpath funciones\
%%
fs = 250;
time = 7.0;
over = 0.1;
TW = 1;
[lags,Twindow] = DFC_timeseg(TW,TW,over,fs,time);
%%
mm = 3
%%

for ss = 1:9
    load(['laplace_BCICIV_2a' filesep 'Lp_folds1_sub' num2str(ss) '.mat'])
    %% Canales 8-12
    X_class11 = Lap{1}{1};
    X_class21 = Lap{1}{2};
    %%
    Aux1 =  size(Lap{1}{1});
    Aux2 =  size(Lap{1}{2});
    lags_Aux = num2cell(lags{1});
    Twindow_Aux = num2cell(repmat(Twindow{1},1,numel(lags_Aux)));
    %
    if   Aux1(1) >= Aux2(1)
        tr_Aux = Aux2(1);
    else
        tr_Aux = Aux1(1);
    end
    X_windClass1FEn = cell(1,tr_Aux);
    X_windClass2FEn = cell(1,tr_Aux);
    %
    X_windClass1PEn = cell(1,tr_Aux);
    X_windClass2PEn = cell(1,tr_Aux);
    %
    X_windClass1SEn = cell(1,tr_Aux);
    X_windClass2SEn = cell(1,tr_Aux);
    
    for tr = 1:tr_Aux
        
        disp(['Sub ' num2str(ss) ' trial ' num2str(tr) ' of 9 .... p1' ])
        X_class1 = cellfun(@(x) squeeze(x) ,mat2cell(X_class11,ones(1,Aux1(1))),'UniformOutput',false);
        X_class2 = cellfun(@(x) squeeze(x) ,mat2cell(X_class21,ones(1,Aux2(1))),'UniformOutput',false);
        X_class1_aux = cellfun(@(p1,p2) X_class1{tr}(:,p1+1:p1+p2),lags_Aux,Twindow_Aux,'UniformOutput',false);
        X_class2_aux = cellfun(@(p1,p2) X_class2{tr}(:,p1+1:p1+p2),lags_Aux,Twindow_Aux,'UniformOutput',false);
        
        %%
        aux1FEn = zeros(22,numel(X_class1_aux));
        aux2FEn = zeros(22,numel(X_class1_aux));
        
        aux1PEn = zeros(22,numel(X_class1_aux));
        aux2PEn = zeros(22,numel(X_class1_aux));
        
        aux1SEn = zeros(22,numel(X_class1_aux));
        aux2SEn = zeros(22,numel(X_class1_aux));
        %
        parfor ii = 1:numel(X_class1_aux)
            %
            x_auxClas1 = num2cell(X_class1_aux{ii},2);
            x_auxClas2 = num2cell(X_class2_aux{ii},2);
            %% Fuzzy Entropy-Mean
            aux1FEn(:,ii) = cell2mat(cellfun(@(x) FuzzyEn2(x,mm,0.75,2),x_auxClas1,'UniformOutput',false));
            aux2FEn(:,ii) = cell2mat(cellfun(@(x) FuzzyEn2(x,mm,0.75,2),x_auxClas2,'UniformOutput',false));
            %% Permutation Entropy
            aux1PEn(:,ii) = cell2mat(cellfun(@(x) PermutationEn(x,mm,1),x_auxClas1,'UniformOutput',false));
            aux2PEn(:,ii) = cell2mat(cellfun(@(x) PermutationEn(x,mm,1),x_auxClas2,'UniformOutput',false));
            %% Sample Entropy
            aux1SEn(:,ii) = cell2mat(cellfun(@(x) SampleEn(x,mm,0.225*std(x)),x_auxClas1,'UniformOutput',false));
            aux2SEn(:,ii) = cell2mat(cellfun(@(x) SampleEn(x,mm,0.225*std(x)),x_auxClas2,'UniformOutput',false));
        end
        X_windClass1FEn{1,tr} = aux1FEn;
        X_windClass2FEn{1,tr} =  aux2FEn;
        
        X_windClass1PEn{1,tr} = aux1PEn;
        X_windClass2PEn{1,tr} = aux2PEn;
        
        X_windClass1SEn{1,tr} = aux1SEn;
        X_windClass2SEn{1,tr} = aux2SEn;
    end
    
    save(['R_Entropy_Channels',filesep,'subj' num2str(ss) '_Tw_1_Over_01.mat'],'X_windClass1FEn','X_windClass2FEn',...
      'X_windClass1PEn','X_windClass2PEn','X_windClass1SEn','X_windClass2SEn')
    %save(['R_Entropy_Channels/subj' num2str(ss) '_Tw_1_Over_09.mat'],'X_windClass1FEn','X_windClass2FEn')
    
end

%% -------------------------------------------------------------------------------------------------------------------------------------
