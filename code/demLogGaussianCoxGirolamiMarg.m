function demLogGaussianCoxGirolamiMarg(rep)

addpath toolbox/; 
dataName = 'logGaussianCoxGirolami_Marg';
storeRes = 1;

for ds=[2 1]
%load the data
% Hyperparameters of model
load ../data/TestData64;
N = 64;
s = 1.91;
b = 1/33;
mu = log(126) - s/2;
m = 1/(N^2);
x_r = 0:1/(N - 1):1;
y_r = 0:1/(N - 1):1;
[xs, ys] = meshgrid(x_r, y_r);
Xt = [];
Xt(:, 1) = xs(:);
Xt(:, 2) = ys(:);

if ds== 2
% work in a subsampled field 
N  = 64/ds;
s  = 1.91;
b  = 1/33;
mu = log(126) - s/2;
m  = 1/(N^2);
x_r = 0:1/(N - 1):1;
y_r = 0:1/(N - 1):1;
[xs, ys] = meshgrid(x_r, y_r);
Xt = [];
Xt(:, 1) = xs(:);
Xt(:, 2) = ys(:);
YY = reshape(Y, 64, 64);
XX = reshape(X, 64, 64);
YY = YY(1:2:64, 1:2:64) + YY(2:2:64, 1:2:64) + YY(1:2:64, 2:2:64) + YY(2:2:64, 2:2:64);
XX = XX(1:2:64, 1:2:64);
Y = YY(:);
X = XX(:);
end

% model options
options = gpsampOptions('LogGaussianCox'); 
options.kern = 'sre';
% create the model
model = gpsampCreate(Y, Xt, options);
model.Likelihood.m = m;
model.Likelihood.logtheta = mu;
model.GP.logtheta = [log((32)*b) log(s)];
%model.GP.logtheta = [0 0];
model.constraints.kernHyper = 'fixed';
model.constraints.likHyper = 'fixed';
model.FF = X(:); 


mcmcoptions.T = 5000;
mcmcoptions.Burnin = 2000;
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
summaryMarg{ds} = summaryStatistics(samples);
summaryMarg{ds}.elapsed = elapsedTime;
summaryMarg{ds}.accRates = accRates;
summaryMarg{ds}.delta = model.delta; 
summaryMarg{ds}.eff_LogL = mcmc_ess(samples.LogL(mcmcoptions.Burnin+1:end));

end

if storeRes == 1
    % don't store samples and model for more than one repeats
    if rep == 1
        save(['../results/' dataName '_repeat' num2str(rep) '.mat'], 'summaryMarg', 'model', 'samples');
    else
        save(['../results/' dataName '_repeat' num2str(rep) '.mat'], 'summaryMarg');
    end
end
