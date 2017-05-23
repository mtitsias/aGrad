function derF = gradlogLIntegralBrownian(lik, Y, F)
%
%Description: Gradient of the Log likelihood wrt mean for the one dimensional 
%             Gaussian distribution   
%

% derivative for standard GP regression
Binter = lik.h*sum(F(1:(end-1)));
sigma2 = exp(lik.logtheta);
derF =  - lik.h*((Binter - lik.mu)/sigma2(1))*ones(length(F),1); 
derF(end) = 0; %- F(end)/0.0001;