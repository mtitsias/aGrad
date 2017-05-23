function [model samples accRates] = gpsampAuxZ_fixedhypers(model, mcmcoptions)
%function [model samples accRates] = gpsampAuxZ_fixedhypers(model, mcmcoptions)
%
% The auxiliary gradient-based method for latent Gaussian models based on the z
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

model.delta = 2/sqrt(n);

if ~isfield(mcmcoptions,'Langevin')
   mcmcoptions.Langevin = 0;
end
Langevin = mcmcoptions.Langevin;

samples.LogL = zeros(1, BurnInIters + Iters);
samples.ff = zeros(1, BurnInIters + Iters);

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:)); 

% If the Langevin algorithm is used then 
% gradient wrt latent variables are required 
if Langevin == 1 
   gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
   derF = gradloglikHandle(model.Likelihood, model.y, F);
   derF2 = derF'*derF;
   meanZ = F + (model.delta/2)*(derF');
end

cnt = 0;

[n D] = size(model.X);
 
twoLambda = 2*model.Lambda;
 
sqrtDeltaLambda = sqrt( (model.delta*model.Lambda)./(twoLambda + model.delta) )'; 
%model.auxPostL = bsxfun(@times, sqrtDeltaLambda',model.U');

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

for it = 1:(BurnInIters + Iters) 
%
    % STEP 1: Gibbs step 
    % sample auxiliary Z
    if Langevin == 1     
        %Z = F + (model.delta/2)*(derF') + sqrt(model.delta/2)*randn(1,n);
        Z = meanZ + sqrt(model.delta/2)*randn(1,n);
    else
        Z = F + sqrt(model.delta/2).*randn(1,n);
    end
     
    % STEP 2: M-H step
    %Fnew = (randn(1, n) + ((2/model.delta)*Z)*model.auxPostL')*model.auxPostL;
    if it <= BurnInIters
        Fnew = ((randn(1, n) + ...
        (((2/model.delta)*Z)*model.U).*sqrtDeltaLambda).*sqrtDeltaLambda)*model.U';
    else
        Fnew = (randn(1, n) + ((2/model.delta)*Z)*model.auxPostL')*model.auxPostL;
    end
    
    % perform an evaluation of the likelihood p(Y | F) 
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));  
    
    % Metropolis-Hastings to accept-reject the proposal
    corrFactor = 0;
    if Langevin == 1 
        
        % new gradient 
        derFnew = gradloglikHandle(model.Likelihood, model.y, Fnew);
        derFnew2 = derFnew'*derFnew;
             
        % correction term
        gzy =  (Z - Fnew)*derFnew - (model.delta/4)*derFnew2;
        gzx =  (Z - F)*derF       - (model.delta/4)*derF2;
        corrFactor = gzy - gzx;  
    end
    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0); 
       
    if accept == 1
         F = Fnew;
         if Langevin == 1
           derF = derFnew;
           derF2 = derFnew2; 
           meanZ = F + (model.delta/2)*(derF');
         end 
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
            sqrtDeltaLambda = sqrt( (model.delta*model.Lambda)./(twoLambda + model.delta) )';
            if Langevin == 1
              meanZ = F + (model.delta/2)*(derF');
            end
        %    
        end
        end
        end
    end
    
    if it == BurnInIters
       model.auxPostL = bsxfun(@times, sqrtDeltaLambda', model.U');
    end
    
    % keep samples after burn in
    if (it > BurnInIters) & (mod(it,StoreEvery) == 0)
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

