function X = FilterWave(x,Wave,lv)
% Esta funcion devuelve los canales filtrados por cada uno de los 2^lv filtros
% Wave: wavelet madre 'sym2';
% lv: nivel de descomposicion  3;

if size(x,1)==1 | size(x,2)==1
    x = x(:);
end

[N N_canal]= size(x);
% Inicializar matriz
X = cell(N_canal,1);
N_fil = 2^lv;
X(:) = {zeros(N,N_fil)};

% filtrar señal
for cnl = 1:N_canal
    % aplicar la desocposiccion de wavelet
    T = wpdec(x(:,cnl),lv,Wave);
    % nodos finales
    nodes = get(T,'tn'); 
    for fil = 1:N_fil
        % reconstruye el nodo fil
        X{cnl}(:,fil) = wprcoef(T,nodes(fil));
    end
end
