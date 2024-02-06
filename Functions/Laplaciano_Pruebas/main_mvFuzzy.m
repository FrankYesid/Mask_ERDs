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
m_s = 22
M   = 3*ones(1,m_s);
r     = 3.0;
n    = 2;
tau = ones(1,m_s);

%%

for ss = 1:9
    
    load(['laplace_BCICIV_2a',filesep,'Lp_folds1_sub' num2str(ss) '.mat'])
    %         y_class1 = ismember(y,1);y_class2 = ismember(y,2);
    %         X_class1 = X(y_class1); X_class2 = X(y_class2);
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
    
    for tr = 1:tr_Aux
        disp(['Sub ' num2str(ss) ' trial ' num2str(tr) ' of 9 .... p1' ])
        X_class1 = cellfun(@(x) squeeze(x) ,mat2cell(X_class11,ones(1,Aux1(1))),'UniformOutput',false);
        X_class2 = cellfun(@(x) squeeze(x) ,mat2cell(X_class21,ones(1,Aux2(1))),'UniformOutput',false);
        X_class1_aux = cellfun(@(p1,p2) X_class1{tr}(:,p1+1:p1+p2),lags_Aux,Twindow_Aux,'UniformOutput',false);
        X_class2_aux = cellfun(@(p1,p2) X_class2{tr}(:,p1+1:p1+p2),lags_Aux,Twindow_Aux,'UniformOutput',false);
        
        %% Fuzzy Entropy
        X_windClass1FEn{1,tr} = cell2mat(cellfun(@(x) mvFE(x,M,r,n,tau),X_class1_aux,'UniformOutput',false));;
        X_windClass2FEn{1,tr} = cell2mat(cellfun(@(x) mvFE(x,M,r,n,tau),X_class2_aux,'UniformOutput',false));;
    end
    
    save(['R_mvFuzzy',filesep,'subj' num2str(ss) '_Tw_1_Over_09.mat'],...
        'X_windClass1FEn','X_windClass2FEn')
    
end

%% -------------------------------------------------------------------------------------------------------------------------------------