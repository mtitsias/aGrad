function out = logLProbit(lik, Y, F)
%
%

yf = Y(:).*F(:); 
     
out = (1+erf(yf/sqrt(2)))/2;         
out = log(out);

