function options = gpsampOptions(gplik) 
%
%

options.kern = 'seard';
options.constraints.kernHyper = 'free'; % 'free' or 'fixed'

switch gplik 
    case 'density'
        % Gaussian likelihood 
        options.Likelihood = 'GPdensity';
    case 'regression'
        % Gaussian likelihood 
        options.Likelihood = 'Gaussian';
    case 'robustRegression'
        % Gaussian likelihood 
        options.Likelihood = 'StudentT';
    case 'classification' % binary classification 
        % Sigmoid likelihood
        options.Likelihood = 'Sigmoid';
    case 'softmax' %  
        % Softmax likelihood
        options.Likelihood = 'Softmax';
    case 'poissonRegr' % Poisson regression
        % Poisson likelihood
        options.Likelihood = 'Poisson';
    case 'offsetPoissonRegr' % Poisson Regression with an offset on the log rate
        % Poisson likelihood after offsetting the log rate
        options.Likelihood = 'OffsetPoisson';
    case 'LogGaussianCox' 
        options.Likelihood = 'LogGaussianCox';  
    case 'IntegralBrownian' 
        options.Likelihood = 'IntegralBrownian'; 
    case 'ODE' 
        % not included 
end

options.constraints.likHyper = 'free'; % 'free' or 'fixed'

        
