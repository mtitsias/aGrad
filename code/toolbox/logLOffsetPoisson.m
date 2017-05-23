function out = logLOffsetPoisson(lik, Y, F)
%function out = logLOffsetPoisson(lik, Y, F)
%
%Description: Poisson GP log-likelihood useful for counts data 
%             This version of the likelihood has an offset for the log rate.
%             Useful if the GP is zero-mean, but we don't believe the process
%             has a zero mean log-rate..
%

offset = lik.logtheta; 
logRate = F(:) + offset;
out = Y(:).*logRate - exp(logRate) - lik.gammaLnY;
