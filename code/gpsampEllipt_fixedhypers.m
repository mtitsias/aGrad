function [model samples] = gpsampEllipt_fixedhypers(model, mcmcoptions)
%function [model samples] = gpsampEllipt_fixedhypers(model, mcmcoptions)
%
% Elliptical slice sampling for latent Gaussian models 
%
% Inputs: 
%         -- model: the structure that contains the log likelihood and the latent
%                    Gaussian prior parameters such as the covariance matrix 
%                    with its pre-computed eigendecomposition (see demos) 
%         -- mcmcoptions: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   state vector of the final MCMC iteration
%                   as well as the learned step size delta
%         -- samples: the structure that contrains the samples 
%
% Michalis K. Titsias (2016)


BurnInIters = mcmcoptions.Burnin; 
Iters = mcmcoptions.T; 
StoreEvery = mcmcoptions.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F = zeros(num_stored, n);
samples.kernLogtheta = zeros(num_stored, model.GP.nParams);
samples.likLogtheta = zeros(num_stored, model.Likelihood.nParams);
samples.LogL = zeros(1, num_stored);
samples.num_calls = 0;
X = model.X;
Y = model.y;
F = model.F; 

samples.LogL = zeros(1, BurnInIters + Iters);
samples.ff = zeros(1, BurnInIters + Iters);

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
%oldL = oldLogLik;
oldLogLik = sum(oldLogLik(:)); 

cnt = 0;
angle_range = 0;

% sample using elliptical slice sampling
for it = 1:(BurnInIters + Iters) 
%
    % the code here follows Iain Murray's implementation   
    nu = randn(1, model.n)*model.L;
    hh = log(rand) + oldLogLik;
    if angle_range <= 0
       % Bracket whole ellipse with both edges at first proposed point
       phi = rand*2*pi;
       phi_min = phi - 2*pi;
       phi_max = phi;
    else
       % Randomly center bracket on current point
       phi_min = -angle_range*rand;
       phi_max = phi_min + angle_range;
       phi = rand*(phi_max - phi_min) + phi_min;
    end

    % Slice sampling loop
    num_calls=0;
    while true
       % Compute xx for proposed angle difference and check if it's on the slice
       Fnew = cos(phi)*F + nu*sin(phi);
       oldLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
       oldLogLik = sum(oldLogLik(:)); 
       num_calls = num_calls + 1;
       if oldLogLik > hh
           % New point is on slice, ** EXIT LOOP **
           break;
       end
       % Shrink slice to rejected point
       if phi > 0
           phi_max = phi;
       elseif phi < 0
           phi_min = phi;
       else
        error('BUG DETECTED: Shrunk to current position and still not acceptable.');
       end
       % Propose new angle difference
       phi = rand*(phi_max - phi_min) + phi_min;
    end
    F = Fnew;
     
       
    samples.num_calls = samples.num_calls + num_calls;
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
    %
        cnt = cnt + 1;
        samples.F(cnt,:) = F;
    %
    end
    %   
    samples.LogL(it) = oldLogLik;
    samples.ff(it) = 0.5*(F*F');    
end
%
%
model.F = F;


