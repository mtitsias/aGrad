function out = logLSoftmax(lik, Y, F)
%function out = logLSoftmax(lik, Y, F)
%
%Description: Softmax GP log-likelihood useful for multiclass classification
%

K = size(F,2);  
M = max(F, [], 2);
out = sum(sum( Y.*F )) - sum(M)  - sum(log(sum(exp(F - repmat(M, 1, K)), 2)));
      
