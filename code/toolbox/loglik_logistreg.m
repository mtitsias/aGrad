function [g, dg, Hessian] = loglik_logistreg(theta, X, Y, T)
%
% 

F = X*(theta(:)); 

YF = -Y.*F;

m = max(0,YF);

g = - sum( m + log( exp(-m) + exp( YF - m )) ); 


if nargout > 1 
%    
  S = sigmoid(F);
  dg = X'*(T - S);
%  
end

if nargout > 2 
%    
    R = S.*(1-S); 
    %Hessian = - X'*(diag(R)*X); 
    Hessian = - X'*bsxfun(@times, R, X);
%    
end


