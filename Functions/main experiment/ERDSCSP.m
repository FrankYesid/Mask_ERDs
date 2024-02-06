%% ERDSCSP

% parametros de la STFT
% Nw = 1; % tama?o de la ventana en seg
Nw = 0.7*fs;
window = hamming(Nw);
noverlap = floor(0.9*Nw);

% Crear matriz de aciertos donde se guardaran todos los resultados por
% fold, por sujeto y por tipo de clasificador
acc = zeros(cv{1}.NumTestSets,numel(X),88*size(X{1}{1},2));

% Tiempo de referencia
t1 = 0;
t2 = 2;

for s=1:numel(1)
    
    %Indices de las etiquetas de interes
    ind = ismember(y{s},labels);
    %Tomar solo las etiquetas de interes
    ys = y{s}(ind);
    fprintf('Subject %d... \n',s)
    tic
    Xs = X{s}(ind);
    
    % Inicializar matrices
    temp = spectrogram(Xs{1}(:,1), window, noverlap,Nw);
    X_suj = zeros(floor((Nw/2)+1),size(temp,2),size(Xs{1},2),numel(Xs));
    
    for tri = 1:numel(Xs)
        for cnl = 1:size(Xs{1},2)
            % signal
            signal = Xs{tri}(:,cnl);
            
            % Calcular STFT
            [temp, f, t] = spectrogram(signal, window, noverlap,Nw,fs);
            
            % Almacenar espectrogramas
            X_suj(:,:,cnl,tri) = abs(temp);
        end
    end
    
    temp1 = abs(t - t1);
    min1 = min(temp1);
    temp2 = abs(t - 2);
    min2 = min(temp2);
    ul = find(temp1 == min1);
    up = find(temp2 == min2);    
    
    for f = 1:cv{s}.NumTestSets
        tr_ind   = cv{s}.training(f); tr_ind = tr_ind(ind);
        ts_ind   = cv{s}.test(f); ts_ind = ts_ind(ind);
        ys_tr = ys(tr_ind);
        X_suj2 = X_suj(:,:,:,tr_ind);
        
        xc1 = X_suj(:,:,:,ys_tr==1);
        xc2 = X_suj(:,:,:,ys_tr==2);
        ERD = zeros(size(xc1,1),size(xc1,2),size(xc1,3),2);
        
        r_nc = squeeze(mean(xc1(:,ul:up,:,:),2));
        r_c = squeeze(mean(r_nc,3));
        m_c = squeeze(mean(xc1,4));
        data = bsxfun(@times,permute(m_c,[2 3 1]),1./r_c) - 1;
        ERD(:,:,:,1) = permute(data,[1 3 2]);
        
        r_nc = squeeze(mean(xc2(:,ul:up,:,:),2));
        r_c = squeeze(mean(r_nc,3));
        m_c = squeeze(mean(xc2,4));
        data = bsxfun(@times,permute(m_c,[2 3 1]),1./r_c) - 1;
        ERD(:,:,:,2) =  permute(data,[1 3 2]);
        
        pm = zeros(size(ERD,2),size(ERD,3));
        
        for ff=1:size(ERD,2)
            for cc = 1:size(ERD,3)
                tem1 = ERD(ff,:,cc,1);
                tem2 = ERD(ff,:,cc,2);
                [~, pm(ff,cc)] = ttest(tem1(:),tem2(:));
            end
        end
        [pvalue, indp] = sort(pm(:),'ascend');
        
        for tri = 1:numel(Xs)
            for cnl = 1:size(Xs{1},2)
                % signal
                signal = Xs{tri}(:,cnl);
                
                % Calcular STFT
                [temp, fa, ti] = spectrogram(signal, window, numel(window)-1,Nw,fs);
                % Almacenar espectrogramas
                X_suja(:,:,cnl,tri) = abs(temp);
            end
        end
        
        ta = 2.5;
        tb = 4.5;
        temp1 = abs(ti - ta);
        min1 = min(temp1);
        temp2 = abs(ti - tb);
        min2 = min(temp2);
        tyt1 = find(temp1 == min1);
        time1 = tyt1(1,1);
        tyt2 = find(temp2 == min2);
        time2 = tyt2(1,1);
        
        Xp = permute(X_suja,[1 3 2 4]);
        Xp = reshape(Xp,size(Xp,1)*size(Xp,2),size(Xp,3),size(Xp,4));
        
        for pp = 6:numel(pvalue)
            Xtmp = Xp(indp(1:pp),:,:);
            C = cell(1,1,numel(Xs));
            parfor trial = 1:numel(Xs)
                % Se computa la covariancia por cada trial.
                C{trial}=cov(Xtmp(:,time1:time2,trial)');
                % Se normaliza la covariancia
                C{trial}=C{trial}/trace(C{trial});
            end
            
            C  = cell2mat(C);
            W        = csp_feats(C(:,:,tr_ind),ys(tr_ind),'train');
            Xc       = csp_feats(C,W,'test');
            Xc_test = Xc(ts_ind,:);
            ys_test = ys(ts_ind);
            mdl = fitcdiscr(Xc(tr_ind,:),ys(tr_ind));
            
            % Clasificar muestras de test usando modelo entrenado.
            acc(f,s,pp) = mean(mdl.predict(Xc(ts_ind,:))==reshape(ys(ts_ind),[sum(ts_ind) 1]));

            if mod(pp,100) == 0
                save ERDSCSP acc 
            end
        end
    end
    
    tictoc = toc;
    fprintf('Acc: %2.1f. Time: %f seconds\n',max(mean(acc(:,s,:),1))*100, tictoc)
    
end
%Guardar los aciertos por sujetos y del clasificador de interés.
save ERDSCSP acc 


