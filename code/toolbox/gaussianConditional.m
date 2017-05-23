function [cmu, cSigma] = gaussianConditional(X, condInds, mu, Sigma)
% Description: It finds the conditional distribution parameters
%              (mean and covariance) of  P(X_rest| X_condInds) 
%
%
% OUTPUT: The conditional mean: cmu, the conditional covariance: cSigma  

d = size(Sigma,1);
d1 = size(condInds(:)',2);
d2 = round(d-d1);

nonCondInds = setdiff([1:d], condInds);
X1 = X(condInds); 
X2 = X(nonCondInds);
mu1 = mu(condInds); 
mu2 = mu(nonCondInds); 

Sigma = Sigma([nonCondInds, condInds],:);
Sigma = Sigma(:,[nonCondInds, condInds]);

B = Sigma(d2+1:end,d2+1:end);
C = Sigma(1:d2, (d2+1):end);

ok = B\C';
cmu = mu2 + (X1 - mu1)*ok;
cSigma = Sigma(1:d2,1:d2) - C*ok;

