function out = logLSigmoid01(lik, Y, F)
%function out = logLSigmoid01(lik, Y, F)
%
%Description: Sigmoid log-likelihood useful for binary classification
%
%           Y is assumed to binary and take values 0-1
%

out = Y(:).*F(:) - log(1 + exp(F(:)));
