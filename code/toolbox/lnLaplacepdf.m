function logpdf = lnLaplacepdf(x,mu,beta)
% natural logarithm of the normal density
%
%

logpdf = -log(2*beta) - (1/beta).*abs(x-mu);