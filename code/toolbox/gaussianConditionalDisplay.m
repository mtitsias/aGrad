function e = gaussianConditionalDisplay(Nc, mu, K)
%
% Illustration of sampling from Gaussian Process conditional 
%

e = 1;

d = size(K,1);


Kb = 0*eye(d);
F = sampleGaussian(1,mu,K+Kb);


M = floor(d/Nc);
obs = [1:M:d]; 
grid = setdiff([1:d], obs);
obs

% compute the conditional Gaussian
[cmu, cSigma] = gaussianConditional(F, obs, mu, K+Kb);

plot(F);
hold on;
plot(obs,F(obs),'ko','MarkerSize',10);

FF = F;
while 1
    FF(grid) = sampleGaussian(1,cmu,cSigma);
    plot(FF,'r')
    plot(F,'b','LineWidth',2);
    plot(obs,F(obs),'go','MarkerSize',10, 'LineWidth', 2);
    pause;
end
