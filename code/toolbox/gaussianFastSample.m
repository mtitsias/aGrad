function X = gaussianFastSample(N, mu, L)
% Function: Generates N samples from a d-dimension
%           Gaussian distribution.
%
% Inputs:  
%       -- N  mumber of samples to generate
%       -- mu 1 x d vector that is the mean of the Gaussian
%       -- L  Choleski decomposition of the covariance matrix (or the square root of 
%          this matrix in case the original covariance was positive semidefinite)
% Outputs: 
%        -- X (N x d) the outputs random vectors  

d = size(L,1);
X = randn(N,d); 
X = X*L;
X = X + ones(N,1)*mu; %repmat(mu,[N,1]);

    