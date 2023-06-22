%% Computing distance between two class of ERD/S
% Key: Adding methods.
% Method 1: inner product.
% Method 2: 1 - inner product.
% Method 3: Euclidean.
% Method 4: Correlation.
% Method 5: 1 - Correlation.
% importante seleccionar un metodo.
% Inputs: 
% erd: Sincronización y desincronización de ambas clases, calculadas anteriormente.
% len:           Tamaño de la ventana analizar.
% over:         Traslape de la ventana traslape si se desea entre.
% nsamples: Numero de muestras.
% nchan:       Numero de canales.
% fs:             Frecuencia de muestro.
% interval:     Intervalo de tiempo para calcular la distancia solicitada.
% sw:            Distancia solicitada.
%                  Distance computation.
% Distance:  1) inner product, 2) 1-inner product, 3) Euclidean,
%                  4) Correlation, 5) 1-Correlation
% 
% Outputs:     Retorna
% Dis_:          La distancia seleccionada por medio de ambas clases.
% Normalizada entre [0 - 1].
% ERDaa:     Los ERD en una .
% Nota: Para el calculo de los ERDs, se analiza el comportamiento en la
% función calcErdsMap() del toolbox de Biosig v3.6
%%
function [Dis_,ERDs] = fncDistcomputation(erd, len,over,nsamples, nchan,fs,interval,sw,distan_,nfreq,MI)
%%
if MI == 0
    t = 1:round(len*fs*over):(nsamples-(len*fs));  % Time.
else
    t = 2.5*fs;
end
labels = size(erd,2);                                           % Numero de la clase.
ERD = cell(1,labels);
ERDs= cell(1,labels);
for cl = 1:labels
    for ch = 1:nchan
        for fre = 1:nfreq
             ERDs{cl}{ch}(fre,:) = erd{cl}.ERDS{ch}.erds(:,fre);
        end
    end
end


for cl = 1:labels
    for ch = 1:nchan
        for fre = 1:nfreq
            for tao = 1:length(t)
                ERD{cl}{ch}{fre,tao} = erd{cl}.ERDS{ch}.erds(t(tao):(t(tao)+len*fs-1),fre);
            end
        end
    end
end

t_ = t./fs;
w = find((t_>interval(1) & t_<interval(2))==1);
nw = numel(w); % time windows in the considered interval
Dis_ = zeros(nchan,nfreq,nw);

if distan_ == 1
    switch sw
        % Calcula la distancia solicitada
        
        case 1
            %% Method 1: inner product
            %         disp('Method 1: inner product')
            %         Dis_ = zeros(nchan,nfreq,nw);
            for ch = 1:nchan
                for fre = 1:nfreq
                    for tao = 1:nw
                        Dis_(ch,fre,tao) = abs(ERD{1}{ch}{fre,w(tao)}'*ERD{2}{ch}{fre,w(tao)}); %norm(ERD_{folds}{1}{ch}{fre}{tao}-ERD_{folds}{2}{ch}{fre}{tao});
                    end
                end
            end
            % Normalization
            % Dis_ = Dis_./max(max(max(Dis_)));
            
        case 2
            %% Method 2: 1 - inner product
            %         disp('Method 2: 1-inner product')
            %         Dis_ = zeros(nchan,nfreq,nw);
            for ch = 1:nchan
                for fre = 1:nfreq
                    for tao = 1:nw
                        Dis_(ch,fre,tao) = 1-abs(ERD{1}{ch}{fre,w(tao)}'*ERD{2}{ch}{fre,w(tao)});
                    end
                end
            end
            % Normalization
            % Dis_ = Dis_./max(max(max(Dis_)));
            
        case 3
            %% Method 3: Euclidean
            %         disp('Method 3: Euclidean')
            %         Dis_ = zeros(nchan,nfreq,nw);
            for ch = 1:nchan
                for fre = 1:nfreq
                    for tao = 1:nw
                        Dis_(ch,fre,tao) = norm(ERD{1}{ch}{fre,w(tao)}-ERD{2}{ch}{fre,w(tao)},2);
                    end
                end
            end
            
            
        case 4
            %% Method 4: Correlation
            %         disp('Method 4: Correlation')
            for ch = 1:nchan
                for fre = 1:nfreq
                    for tao = 1:nw
                        Dis_(ch,fre,tao) = abs(corr(ERD{1}{ch}{fre,w(tao)},ERD{2}{ch}{fre,w(tao)}));
                    end
                end
            end
            
        case 5 
            %% Method 5: 1 - Correlation
            %         disp('Method 5: Correlation')
            for ch = 1:nchan
                for fre = 1:nfreq
                    for tao = 1:nw
                        Dis_(ch,fre,tao) = 1-abs(corr(ERD{1}{ch}{fre,w(tao)},ERD{2}{ch}{fre,w(tao)}));
                    end
                end
            end
            
        otherwise
            disp('Distance 1) inner product, 2) 1-inner product, 3) Euclidean, 4) Correlation, 5) 1-Correlation')
            
    end
    
    % Normalization
    Dis_ = Dis_./max(max(max(Dis_)));
end
%%Graphs
%imagesc(squeeze(Dis_corr(:,3,:)))
