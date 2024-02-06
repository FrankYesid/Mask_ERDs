function FEn = FuzzyEn_mean(X,m,r)
N = length(X);
phi = zeros(1,2);
%% Normalization -- Step 1
 X_n = zscore(X);
 %X_n = X;
%%
 %patterns = zeros(m+1,N-m);
% for i = 1:m+1
%     patterns(i,:) = X_n(i:N-m+i-1);
% end
%% Step 2
for mm = m:m+1
    %% Step 3
    patterns = zeros(mm,N-mm+1);
    for i = 1:mm
        patterns(i,:) = X_n(i:N-mm+i);
    end
    patterns = patterns - mean(patterns);
    dist = pdist(patterns','chebychev'); 
    Dg = exp(-log(2)*((dist/r).^2));
    %Dg = (-((dist.^2)/r));
    %% Step 4 and 5
    phi(mm-m+1) = (sum(Dg))/(N-mm-1)/(N-mm);
end
%% Fuzzy Entropy
FEn = log(phi(1)) -log(phi(2));
end