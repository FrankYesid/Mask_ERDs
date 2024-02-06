function [Ktr,Kts,sigma]=gaussKernel(Xtr,Xts)

Dtr     = pdist(Xtr);
sigma   = median(Dtr);
Ktr     = exp(-squareform(Dtr).^2/(2*sigma^2));
Kts     = exp(-pdist2(Xtr,Xts).^2/(2*sigma^2));
