function model = gpsampCreate(y, X, options, flag)
%
%

if nargin == 3
    flag = 0;
end

model.type = 'gpmodel';
model.Likelihood.type = options.Likelihood;

[numData D] = size(X);

% if you have multivariate outputs 
J = size(y,2); 

model.y = y; 
model.X = X;
model.numData = numData;
model.D = D;
model.n = numData;
model.J = J;
switch model.Likelihood.type 
    case 'GPdensity' % standard regression
         %
         model.Likelihood.nParams = 1; % parameters (excluding gp function F)
         model.Likelihood.p0.mu = mean(X); 
         model.Likelihood.p0.sigma2 = 10*max(var(X));  
         model.Likelihood.X = X;
    case 'Gaussian' % standard regression
         %
         model.Likelihood.nParams = 1; % parameters (excluding gp function F)
         model.Likelihood.logtheta = log(0.05); % log(sigma2) 
    case 'StudentT' % robust regression
         %
         % the degrees of freedom are treated as fixed to value 4
         model.Likelihood.v = 4;
         model.Likelihood.nParams = 1;  
         % variance parameter
         model.Likelihood.logtheta = log(0.05); % log(sigma2) 
         %
    case 'Probit'  % binary classifcation
         %
         model.Likelihood.nParams = 0;
         %
    case 'Sigmoid' % binary classification
         %
         model.Likelihood.nParams = 0;
         %
    case 'Softmax' % multiclass classification
         %
         model.Likelihood.nParams = 0;
         %   
    case 'Poisson' % for counts data      
         %
         model.Likelihood.nParams = 1; 
         model.Likelihood.logtheta = 0;
         model.Likelihood.gammaLnY = gammaln(y(:)+1);
         %
    case 'OffsetPoisson' % for counts data with offset on log rate
         %
         model.Likelihood.nParams = 1;
         model.Likelihood.logtheta = 0;
         model.Likelihood.gammaLnY = gammaln(y(:)+1);
         
    case 'LogGaussianCox' % log-Gaussian Cox likelihood in regular grid
         %
         model.Likelihood.nParams = 1;
         model.Likelihood.logtheta = 0;
         % area of each small patch in the discretized grid 
         model.Likelihood.m = 1; 
    case 'IntegralBrownian' % log-Gaussian Cox likelihood in regular grid
         %
         model.Likelihood.nParams = 1;
         model.Likelihood.logtheta = 0;
         model.Likelihood.h = 0.1;
         %
    case 'ODE'
         % %%
end     

% sample the hyperparameters or keep then fixed
model.constraints.kernHyper = options.constraints.kernHyper;
model.constraints.likHyper = options.constraints.likHyper;
%model.singleScaleParam = 0;
%model.alpha = 1;

if (J == 1) && (flag == 0) 
%     
model.GP.jitter = 0;
switch options.kern 
    case 'se'
       model.GP.type = 'se';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP.logtheta = [2*mean(log((max(X) - min(X))*0.1)) 0];
       model.GP.nParams = 2;
    case 'seard'
       model.GP.type = 'seard';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP.logtheta = [2*log((max(X) - min(X))*0.1) 0];
       model.GP.nParams = D+1;
    case 'sre'
       model.GP.type = 'sre';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP.logtheta = [2*mean(log((max(X) - min(X))*0.1)) 0];
       model.GP.nParams = 2;
    case 'sreard'    
       model.GP.type = 'sreard';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP.logtheta = [2*log((max(X) - min(X))*0.1) 0];
       model.GP.nParams = D+1;   
end
%
else
%    
  for j=1:model.J
  model.GP{j}.jitter = 0;
  switch options.kern 
    case 'se'
       model.GP{j}.type = 'se';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP{j}.logtheta = [2*mean(log((max(X) - min(X))*0.1)) 0];
       model.GP{j}.nParams = 2;
    case 'seard'
       model.GP{j}.type = 'seard';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP{j}.logtheta = [2*log((max(X) - min(X))*0.1) 0];
       model.GP{j}.nParams = D+1;
    case 'sre'
       model.GP{j}.type = 'sre';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP{j}.logtheta = [2*mean(log((max(X) - min(X))*0.1)) 0];
       model.GP{j}.nParams = 2;
    case 'sreard'    
       model.GP{j}.type = 'sreard';
       % kernel hyperparameters (lengthscales, kernel variance) in the log space
       model.GP{j}.logtheta = [2*log((max(X) - min(X))*0.1) 0];
       model.GP{j}.nParams = D+1;   
  end  
  end
%    
end

% prior over the likelihood parameters  
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.likParams.type = 'normal';
model.prior.likParams.constraint = 'positive';
model.prior.likParams.priorSpace = 'log';
model.prior.likParams.a = 0; % mean 
model.prior.likParams.b = 1; % variance
 
model.jitter = 1e-6;

% prior over GP kernel hyperparameters 
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.kernParams.type = 'normal';
model.prior.kernParams.constraint = 'positive';
model.prior.kernParams.priorSpace = 'log';
model.prior.kernParams.a = 0; % mean 
model.prior.kernParams.b = 1; % variance
 
% GP latent function values needed to define the likelihood
model.F = zeros(model.J,model.numData);
