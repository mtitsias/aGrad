function out = lnnormalpdf(x,mu,sigma2, grad)
% natural logarithm of the normal density
%
%

if nargin == 3
    grad = 0;
end

if grad == 0
  out = -0.5*log(2*pi*sigma2) - (0.5/sigma2).*((x-mu).^2);
else
  % return also the gradient wrt x
  out = (mu-x)/sigma2;
  out = out(:); 
end