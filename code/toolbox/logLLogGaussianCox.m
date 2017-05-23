function out = logLLogGaussianCox(lik, Y, F)
%function out = logLLogGaussianCox(lik, Y, F)
%
%Description: Poisson GP log-likelihood useful for counts data 
%

offset = lik.logtheta(1); 
logRate = F(:) + offset;
out = Y(:).*logRate - lik.m*exp(logRate);