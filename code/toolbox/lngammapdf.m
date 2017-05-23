function logpdf = lngammapdf(x,a,b)
%
%

%z = x ./ b;
%u = (a - 1) .* log(z) - z - gammaln(a);
%y = exp(u) ./ b;

logpdf = (a-1)*log(x) - b*x + a*log(b) - gammaln(a);