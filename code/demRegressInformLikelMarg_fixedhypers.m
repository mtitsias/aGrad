function demRegressInformLikelMarg_fixedhypers(rep)

addpath toolbox/;

dataName = 'regressionInformLikeToy_Marg';
storeRes = 1;

for i=1:3

% Fix seeds
%randn('seed', 1e5);
%rand('seed', 1e5);
%sigma2vals = [1^2 0.1^2 0.01^2 0.001^2];
%step = 0.004; 
%X = (0:step:4)';
%sigma2 = sigma2vals(i); % 0.1^2; 
%Y = sin(2*pi*X) + cos((2.5*pi*X)) + sqrt(sigma2)*randn(size(X,1),1);
%XX{i} = X;
%YY{i} = Y;
%save data/regressinformlik.mat XX YY sigma2vals;
load ../data/regressinformlik.mat;
X = XX{i};
Y = YY{i};
sigma2 = sigma2vals(i);

% model options
options = gpsampOptions('regression'); 

% create the model
model = gpsampCreate(Y, X, options);
model.Likelihood.logtheta = log(sigma2);

model.constraints.kernHyper = 'fixed';
model.constraints.likHyper = 'fixed';

mcmcoptions.T = 5000;
mcmcoptions.Burnin = 10000;
mcmcoptions.StoreEvery = 1;
mcmcoptions.Langevin = 1;

% precompute the inverse covariacne matrix Q 
model.K = kernCompute(model.GP, model.X);
[model.U, model.Lambda, tmp] = svd(model.K);
model.Lambda = diag(model.Lambda);

tic;
[model samples accRates] = gpsampAuxMarg_fixedhypers(model, mcmcoptions);
elapsedTime = toc;

% compute statistics 
summaryMarg{i} = summaryStatistics(samples);
summaryMarg{i}.elapsed = elapsedTime;
summaryMarg{i}.accRates = accRates;
summaryMarg{i}.delta = model.delta; 
summaryMarg{i}.eff_LogL = mcmc_ess(samples.LogL(mcmcoptions.Burnin+1:end));

end

if storeRes == 1
    save(['../results/' dataName '_repeat' num2str(rep) '.mat'], 'summaryMarg');
end

