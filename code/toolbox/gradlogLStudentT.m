function derF = gradlogLStudentT(lik, Y, F)
%function derF = gradlogLStudentT(lik, Y, F)
%
%

v = lik.v; 
sigma2 = exp(lik.logtheta(1)); 

diff = Y(:)-F(:);
num = 1 + (1/(v*sigma2))*(diff.^2);

derF = ((v+1)/(v*sigma2))*(diff./num);
 