
%clear all; 
close all; 
addpath toolbox/;
addpath ../data/;
 
dataName = 'mnistsoftmax';

% Load the MNIST dataset and 
% create the appropritate input and output data matrices 
load mnist_all.mat;
% number of classes
K = 10;
T = []; 
X = [];
for j=1:10
% 
    s = ['train' num2str(j-1)];
    Xtmp = eval(s); 
    Xtmp = double(Xtmp(1:100,:));   
    Ntrain = size(Xtmp,1);
    Ttmp = zeros(Ntrain, K); 
    Ttmp(:,j) = 1; 
    X = [X; Xtmp]; 
    T = [T; Ttmp]; 
%    
end

% normalize the pixels to take values in [0,1]
X = X/255;     
Y = T;
   
[n D] = size(X);
   
% model options
options = gpsampOptions('softmax'); 
options.kern = 'se';
   
model = gpsampCreate(Y, X, options, 1); 
for j=1:model.J 
    model.GP{j}.jitter = 1e-8;
    model.GP{j}.logtheta = [2*max(log((max(X) - min(X))*0.5*sqrt(D))) 0];
end
model.XX = -2*(X*X') + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   
model.constraints.kernHyper = 'free';
model.constraints.likHyper = 'free';
   
mcmcoptions.T = 5000;
mcmcoptions.Burnin = 50000;
mcmcoptions.StoreEvery = 1;
mcmcoptions.Langevin = 1;

tic;
[model samples accRates] = gpsamppCN_learnhypers(model, mcmcoptions);
elapsedTime = toc;

% compute statistics 
eff_F = []; 
eff_kern = []; 
for j=1:model.J
   samp = samples;     
   samp.F = squeeze(samples.F(j,:,:))';
   tmp = summaryStatistics(samp);
   eff_F = [eff_F, tmp.eff_F];
   tmp = [];
   for k=1:model.GP{j}.nParams 
      tmp(k) = mcmc_ess(samples.kernLogtheta{j}(:,k));
   end
   eff_kern = [eff_kern, tmp];
end
summarypCNL.elapsed = elapsedTime;
summarypCNL.accRates = accRates;
summarypCNL.delta = model.delta; 
summarypCNL.kerndelta = model.kerndelta;
summarypCNL.eff_F = eff_F;
summarypCNL.eff_kern = eff_kern;
summarypCNL.LogL = samples.LogL;
 
save(['../results/' dataName '_pCNL_learnhypers.mat'], 'summarypCNL', 'model', 'samples');

