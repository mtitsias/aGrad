function out = logLStudentT(lik, Y, F)
%function out = logLStudentT(lik, Y, F)
%
%Description: Log likelihood for the one dimensional 
%             Gaussian distribution   
%

%v = exp(lik.logtheta(1)); 
v = lik.v; 
sigma2 = exp(lik.logtheta(1)); 

out = - 0.5*log(2*pi*v*sigma2) + gammaln(0.5*(v+1)) - gammaln(0.5*v); 

out = out - 0.5*(v+1)*log(1  + (1/(v*sigma2))*((Y(:)-F(:)).^2) ); 

