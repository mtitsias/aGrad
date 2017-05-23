function [psi0, Psi2] = kernVardistPsiCompute(model)
%
%


[n D] = size(model.X);

p0 = model.Likelihood.p0;

Z = model.X;

% psi0
psi0 = 1;


alpha = 1./exp(model.GP.logtheta(1));

AS_n = (1 + 2*alpha*p0.sigma2)^0.5;  
   
%Z_n = (repmat(vardist.means(n,:),[M 1]) - Z)*0.5;
Z_n = bsxfun(@minus, p0.mu, Z)*0.5;
%Z_n = Z_n.*repmat(sqrt(A)./AS_n,[M 1]);
Z_n = bsxfun(@times, Z_n, sqrt(alpha)./AS_n);
distZ = dist2(Z_n,-Z_n);    
% ZZ = Z.*(repmat(sqrt(A),[M 1]));
ZZ = bsxfun(@times, Z, sqrt(alpha));
distZZ = dist2(ZZ,ZZ);

% dont multiply with sigmaf^4 
Psi2 = exp(-0.25*distZZ - distZ  - 0.5*D*log(1 + 2*alpha*p0.sigma2));

