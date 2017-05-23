function [model samples accRates] = gpsampAuxU_fixedhypers(model, mcmcoptions)
%function [model samples accRates] = gpsampAuxU_fixedhypers(model, mcmcoptions)
%
% The auxiliary gradient-based method for latent Gaussian models based on the u
% auxiliary variable 
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
Y = model.y;
F = model.F; 


samples.LogL = zeros(1, BurnInIters + Iters);
samples.ff = zeros(1, BurnInIters + Iters);

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:)); 

model.delta = 2/sqrt(n); 

twoLambda = 2*model.Lambda;
Lambda1 = ( (model.delta*model.Lambda)./(twoLambda + model.delta) )';
sqrtDeltaLambda = sqrt(Lambda1);

% gradient wrt latent variables are required 
gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
derF = gradloglikHandle(model.Likelihood, model.y, F);
derFU = model.U'*derF;

cnt = 0;

[n D] = size(model.X);
         
range = 0.05; 
opt = 0.54;

% adaption step size
epsilon = 0.05;

acceptHistF = zeros(1, BurnInIters + Iters);

for it = 1:(BurnInIters + Iters) 
%
    % STEP 1: Gibbs step 
    % sample auxiliary U
    U = F + sqrt(model.delta/2).*randn(1,n);
     
    twodeltaU = (2/model.delta)*U; 
    
    % STEP 2: M-H step
    % sample new state vector Fnew given old F (y given x in paper's notation)
    UtwodeltaU = twodeltaU*model.U; 
    Fnew = (randn(1, n).*sqrtDeltaLambda + (derFU' + UtwodeltaU).*Lambda1)*model.U';
   
    % perform an evaluation of the likelihood p(Y | F) 
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));  
    
    % Metropolis-Hastings to accept-reject the proposal 
    % new gradient 
    derFnew = gradloglikHandle(model.Likelihood, model.y, Fnew);
    derFUnew = model.U'*derFnew;
    rhoxyu = F*derFnew  -  ( Lambda1.*(UtwodeltaU + 0.5*derFUnew'))*derFUnew; 
    rhoyxu = Fnew*derF  -  ( Lambda1.*(UtwodeltaU + 0.5*derFU'))*derFU;     
    corrFactor = rhoxyu - rhoyxu; 
    
    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0); 
          
    if accept == 1
         F = Fnew;
         derF = derFnew;
         derFU = derFUnew; 
         oldLogLik = newLogLik; 
    end
    acceptHistF(it) = accept;    
    
    % Adapt proposal during burnin 
    if mod(it,5) == 0 
        if (it >= 50) 
        accRateF = mean(acceptHistF((it-49):it))*100;
        if (it <= BurnInIters)
        if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
        %    
            model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;  
            Lambda1 = ( (model.delta*model.Lambda)./(twoLambda + model.delta) )';
            sqrtDeltaLambda = sqrt(Lambda1);
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
    samples.LogL(it) = oldLogLik;       
    samples.ff(it) = 0.5*sum(sum(F.*F));   
    %        
end
%
%

model.F = F;
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100; 
