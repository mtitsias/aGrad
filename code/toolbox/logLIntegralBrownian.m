function out = logLIntegralBrownian(lik, Y, F)
%
%Description: Log likelihood for the one dimensional 
%             Gaussian distribution   
%

sigma2 = exp(lik.logtheta); 


% there is only one observation 
Binter = lik.h*sum(F(1:(end-1)));
out = - 0.5*log(2*pi*sigma2(1)) - (0.5/sigma2(1))*((Binter - lik.mu)^2);
%out = out - 0.5*log(2*pi*0.0001) - (0.5/0.0001)*(F(end)^2);
