function derF = gradlogLLogGaussianCox(lik, Y, F)
%function derF = gradlogLLogGaussianCox(lik, Y, F)
%
%Description: Derivative of the Poisson GP log-likelihood useful for counts data 
%


derF = Y(:) - (lik.m*exp(lik.logtheta(1)))*exp(F(:));         