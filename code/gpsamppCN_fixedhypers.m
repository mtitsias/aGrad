function [model samples accRates] = gpsamppCN_fixedhypers(model, mcmcoptions)
%function [model samples accRates] = gpsamppCN_fixedhypers(model, mcmcoptions)
%
% pCP for latent Gaussian models  
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
%         -- accRates: acceptance rate
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

if ~isfield(mcmcoptions,'Langevin')
   mcmcoptions.Langevin = 0;
end
Langevin = mcmcoptions.Langevin;

samples.LogL = zeros(1, BurnInIters + Iters);
samples.ff = zeros(1, BurnInIters + Iters);

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
%oldL = oldLogLik;
oldLogLik = sum(oldLogLik(:)); 

model.delta = 2/sqrt(n);
model.beta = 2/(model.delta + 2);
beta = model.beta;
sqrtOnebeta2 = sqrt(1-beta^2);

% If the Langevin algorithm is used then 
% gradient wrt latent variables are required 
if Langevin == 1 
   gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
   derF = gradloglikHandle(model.Likelihood, model.y, F);
   CderF = model.K*derF;
end

cnt = 0;

range = 0.05;
if isfield(model,'opt')
   opt = model.opt;
   range = 0; 
else
   if Langevin == 1 
      opt = 0.54;
   else
      opt = 0.25;
   end
end
% adaption step size
epsilon = 0.05;

acceptHistF = zeros(1, BurnInIters + Iters);

% sample using Neal's method
for it = 1:(BurnInIters + Iters) 
%     
    if Langevin == 1     
        Fnew = beta*F + (1-beta)*(CderF') + sqrtOnebeta2*(randn(1, model.n)*model.L);
    else
        Fnew = beta*F + sqrtOnebeta2*(randn(1, model.n)*model.L);
    end
    
    % perform an evaluation of the likelihood p(Y | F) 
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));  
     
    % Metropolis-Hastings to accept-reject the proposal
    if Langevin == 1 
        
        % new gradient 
        derFnew = gradloglikHandle(model.Likelihood, model.y, Fnew);
        CderFnew = model.K*derFnew;
         
        hxy = newLogLik + (1/(1+beta))*( (F - beta*Fnew - (0.5*(1-beta))*CderFnew')*derFnew );
        hyx = oldLogLik + (1/(1+beta))*( (Fnew - beta*F - (0.5*(1-beta))*CderF')*derF );
    else
        hxy = newLogLik;
        hyx = oldLogLik;
    end
     
    [accept, uprob] = metropolisHastings(hxy, hyx, 0, 0); 
   
    acceptHistF(it) = accept;  
     
    if accept == 1
        F = Fnew;   
        if Langevin == 1
          derF = derFnew;
          CderF = CderFnew;
        end
        oldLogLik = newLogLik; 
    end
    
    % Adapt proposal during burnin
    if mod(it,5) == 0 
        if (it >= 50) 
        accRateF = mean(acceptHistF((it-49):it))*100;
        if (it <= BurnInIters)
        if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
        %    
            model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;
            model.beta = 2/(model.delta + 2);
            beta = model.beta;
            sqrtOnebeta2 = sqrt(1-beta^2);
        %    
        end
        end
        end
    end
    
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
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100; 

 
