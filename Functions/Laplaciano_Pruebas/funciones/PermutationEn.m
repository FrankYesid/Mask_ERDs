function PEn = PermutationEn(X,Order,Delay)
% Inputs
% X  = señal de entrada X[nx1]; n : tamaño de la señal
% Order = embedding dimension controla la longuitud de  cada nuevo vector
% fila (entero)
% Delay = embedding time delay controla la numero de periodos de tiempo
% entre lo elemnentos de cada uno de los nuevo vectores fila (entero)
% Alfa  =  [-1<= Alfa <=1], Alfa = 0 se obtiene la Permutation Entropy (vector de float)

m = Order;
tau = Delay;

%% Permutation Entropy convesinal
n = length(X);
K = n-(m-1)*tau;
Xk = zeros(K,m);

for ii = 1:K
    Xk(ii,:) = X(ii:tau:ii+(m-1)*tau);
end

[~,sort_indMat] =sort(Xk,2);
indMat = zeros(size(sort_indMat));
for ii = 1:size(sort_indMat,1)
    indMat(ii,sort_indMat(ii,:)) = 1:m;
end

indVec = m.^(0:m-1) * (indMat'-1);

[~,ia,~] = unique(sort(indVec),'first');
nSymbol = diff([ia;length(indVec)+1]);

pSymbol = nSymbol/sum(nSymbol);
PEn = -sum(pSymbol.*log2(pSymbol)); %PE

end

