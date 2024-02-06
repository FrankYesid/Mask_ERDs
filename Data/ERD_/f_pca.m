% % Analisis de relevancia
% clear all; clc;
% load('X_sujetos_Izq_Der.mat')

function Sujetos = f_pca(X_sujetos_Izq_Der)
%X_sujetos_Izq_Der{sujetos}(trial,canal,frecuencia,tiempo)

Sujetos = cell(9,1);
for su = 1:9;
    fprintf(['sujeto: ' num2str(su) ' de ' '9' '\n'])
    [N_tri N_canal N_freq N_time] = size(X_sujetos_Izq_Der{su});
    Z = cell(22,1);
    S = cell(22,1);
    Z(:) = {zeros(N_tri, N_freq*N_time)};
    S(:) = {zeros(N_freq,N_time)};
    
    for cnl = 1:N_canal
        for tri = 1:N_tri;
            temp = squeeze(X_sujetos_Izq_Der{su}(tri,cnl,:,:));
            % vectorizar tiempo frecuencia
            Z{cnl}(tri,:) = temp(:)';
        end
        % calcular la magnitud y estandarizar
        Z{cnl} = zscore(abs(Z{cnl}));
        %PCA
        [w, ~, lambda, tsquared, explained] = pca(Z{cnl},'NumComponents',N_tri);
        rho_f = abs(w)*lambda(1:size(w,2));
        % reconstruir la matriz
        S{cnl} = reshape(rho_f,N_freq,N_time);
    end
    Sujetos{su} = S;
end
