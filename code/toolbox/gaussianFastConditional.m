function [cmu, cSigma, KInvK] =  gaussianFastConditional(mu, K, Express, Given)
% Partial computation of the conditional X(Express) | X(Given)
% Inputs: 
%      -- mu the mean of the Gaussian as a row vector 
%      -- K the covariance matrix of the Gaussian  
%      -- Express and Given are disjoint sets of indices pointing 
%         in the random vector 
%  
% The function writes the Gaussian N(mu,K) in the form 
%
%   [X(Express)]      mu(Express)  [ A  C  ...]
%   [X(Given)  ]  = N(mu(Given),   [ C^T B ...])
%   [X(Rest)   ]      mu(Rest)     [ ...   ...]
% 
% and then computes elements of the conditional:
%        N(mu(Express) + C*B^(-1)*[X(Given) - mu(Given)], A - C*B^(-1)*C') 
%
% Outputs:
%       -- cmu    =  mu(Express)' - mu(Given)'*(B^(-1)*C') 
%       -- cSigma =  A - C*B^(-1)*C'
%       -- KInvK  =  B^(-1)*C'
%

jitter = 10e-8;

SizF = size(mu,2);
SizE = size(Express,2);
SizG = size(Given,2);
Ignore = setdiff([1:SizF], [Express, Given]);


muE = mu(Express); 
muG = mu(Given);

K = K([Express, Given, Ignore],:);
K = K(:,[Express, Given, Ignore]);
K = K(1:SizE+SizG,1:SizE+SizG);

B = K(SizE+1:end,SizE+1:end) + jitter*K(1,1)*eye(SizG);
C = K(1:SizE, (SizE+1):end);

KInvK = B\C';
cmu = muE - muG*KInvK;
cSigma = K(1:SizE,1:SizE) - C*KInvK;
