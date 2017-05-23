function out = logLPoisson(lik, Y, F)
%function out = logLPoisson(lik, Y, F)
%
%Description: Poisson GP log-likelihood useful for counts data 
%

offset = lik.logtheta; 
logRate = F(:) + offset;
out = Y(:).*logRate - exp(logRate) - lik.gammaLnY;

