function derF = gradlogLSigmoid(lik, Y, F)
%function derF = gradlogLSigmoid(lik, Y, F)
%
%Description: Derivative of sigmoid GP log-likelihood useful for binary classification
%

% derivative for binary classification  
derF = 0.5*(Y(:)+1) - 1./(1+exp(-F(:)));  % derivative for binary classification  