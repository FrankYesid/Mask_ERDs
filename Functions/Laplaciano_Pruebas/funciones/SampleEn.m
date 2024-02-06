function  SEn  = SampleEn( X,m,r )
% X -time serie [1 x N]
% m -embedding dimension (int)
% r -the distance between x(i) and x(j) lower than the tolerance
%     r = (0.1 ~ 0.25)*std(X)
%% Sample Entropy
% -------------------State Space Reconstruction
N = numel(X);
m_aux = m;
i = 1:N -m_aux +1;
x_m = zeros(numel(i),m_aux);
for ii = i
    x_m(ii,:) = X(ii:ii+m_aux-1);
end
% --------------------Distance
Dist = pdist(x_m, 'chebychev');
%Dist = pdist(x_m, 'euclidean');
% ---------------Calculate Bm
Bm1 = 2*sum(Dist<=r)/(N-m_aux-1)/(N-m_aux);
%Bm1 = 2*sum(Dist<=r)/(N-m_aux+1)/(N-m_aux);
%%
% ---------------------------------------------------------
clearvars x_m Dist i
% ---------------------------------------------------------
%%
m_aux = m+1;
i = 1:N-m_aux +1;
x_m = zeros(numel(i),m_aux);
for ii = i
    x_m(ii,:) = X(ii:ii+m_aux-1);
end
% --------------------Distance
Dist = pdist(x_m, 'chebychev');
%Dist = pdist(x_m, 'euclidean');
% ---------------Calculate Bm
Bm2 = 2*sum(Dist<=r)/(N-m_aux-1)/(N-m_aux);
%Bm2 = 2*sum(Dist<=r)/(N-m_aux+1)/(N-m_aux);
%%
% -------- sample Entropy
SEn = -log(Bm2/Bm1);
%%
% If A=0 or B=0, SampEn would return an infinite value. However, the
% lowest non-zero conditional probability that SampEn should
% report is A/B = 2/[(N-m-1)(N-m)]
if isinf(SEn) || isnan(SEn)
    % Note: SampEn has the following limits:
    %       - Lower bound: 0
    %       - Upper bound: log(N-m)+log(N-m-1)-log(2)
    SEn = -log(2/((N-m-1)*(N-m)));
end
end

