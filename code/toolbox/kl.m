function out = KL(mu0, Sigma0, mu1, Sigma1)
% function out = KL(mu0, Sigma0, mu1, Sigma1)
%
% Description :  Computes the  KL divergence between two Gaussians distribution 
%                
%

mu0 = mu0(:);
mu1 = mu1(:);

N = size(Sigma0,1);

L0 = jitterChol(Sigma0)'; 
L1 = jitterChol(Sigma1)'; 
invSigma1 = L1'\(L1\eye(N));

out = sum(log(diag(L1))) - sum(log(diag(L0)))...
      + 0.5*trace(invSigma1*Sigma0) + 0.5*((mu1-mu0)'*(invSigma1*(mu1-mu0)))...
      - 0.5*N;
