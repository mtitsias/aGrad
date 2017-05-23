function out = logLGaussian(lik, Y, F)
%function out = logLGaussian(Y,F,sigma2)
%
%Description: Log likelihood for the one dimensional 
%             Gaussian distribution   
%

sigma2 = exp(lik.logtheta); 

n = size(Y(:),1);
out = -0.5*n*log(2*pi*sigma2) - (0.5/sigma2)*sum((Y(:)-F(:)).^2);
