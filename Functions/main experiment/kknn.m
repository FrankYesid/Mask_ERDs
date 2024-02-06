function acc=kknn(Xtr,ytr,Xts,yts,nn,input)

if strcmp(input,'features')    
    Dtr     = pdist(Xtr);
    sigma   = median(Dtr);
    Kts     = exp(-pdist2(Xtr,Xts).^2/(2*sigma^2));
    
elseif strcmp(input,'kernels')
    Kts = Xts;
end

[~,ind] = sort(Kts,'descend');
acc     = mean(mode(ytr(ind(1:nn,:)),1)' == yts);
