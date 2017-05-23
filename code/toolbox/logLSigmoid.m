function out = logLSigmoid(lik, Y, F)
%function out = logLSigmoid(lik, Y, F)
%
%Description: Sigmoid GP log-likelihood useful for binary classification
%

YF= Y(:).*F(:); 
out = -log(1 + exp(-YF));
