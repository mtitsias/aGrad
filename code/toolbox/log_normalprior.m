function [g, dg, Hessian] = log_normalprior(Theta, m, Sigma)
%




D = length(Theta);  

d = Theta(:) - m;  
   
g = - 0.5*D*log(2*pi) - 0.5*sum(log(Sigma)) - 0.5*sum((d.^2)./Sigma);
   
if nargout > 1 
   dg = - d./Sigma;  
end
   
if nargout > 2 
   Hessian = - diag(1./Sigma);  
end