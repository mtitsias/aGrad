function out = logLGPdensity(lik, y, F)


X = lik.X;

[n D] = size(X);

mu = lik.p0.mu; 
sigma2 = lik.p0.sigma2;
D = length(mu); 
logp0 = - 0.5*D*n*log(2*pi*sigma2) - (0.5/sigma2)*sum(sum( ( X -  repmat(mu,n,1) ).^2));

out = logp0  + sum( log( (F.^2)) );

FinvL = (F(:)')*lik.invL;
out = out - n*log(FinvL*(lik.C*FinvL')   + lik.psi0PlusTrace );
