function derF = gradlogLSoftmax(lik, Y, F)
%function derF = gradlogLSoftmax(lik, Y, F)
%
%Description: Derivative of softmax GP log-likelihood useful for multiclass classification
%

% derivative for softmax classification  
derF = Y - softmax(F); 