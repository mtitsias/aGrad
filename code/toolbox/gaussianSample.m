function X = gaussianSample(N, mu, Sigma)
%Description:  It generates N samples from a d-dimensional
%              Gaussian distribution.
%
%N  : Number of samples to generate
%mu : mean
%Sigma : covariance matrix for the samples, should be positive definite

d = size(Sigma,1);

X=randn(N,d); 

% Cholesky decomposition for a positive definite Covariance
% If the covariance is not positive definite uses the square root 
[L,r]=chol(Sigma);
if r>0 
    X = real(X*sqrtm(Sigma));
else
    X = X*L;
end
X=X+repmat(mu,[N,1]);

    