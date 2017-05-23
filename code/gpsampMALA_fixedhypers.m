function [model samples accRates] = gpsampMALA_fixedhypers(model, mcmcoptions)
%function [model samples accRates] = gpsampMALA_fixedhypers(model, mcmcoptions)
%
% pMALA for latent Gaussian models  
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
samples.LogL = zeros(1, num_stored);
samples.Log_joint = zeros(1, num_stored);
samples.ff = zeros(1, num_stored);
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
oldLogLik = sum(oldLogLik(:)); 
QF = F*model.Q;
QF2 = QF*F';
oldLog_joint = oldLogLik - 0.5*QF2;

model.delta = 1/sqrt(n); 

% If the Langevin algorithm is used then 
% gradient wrt latent variables are required 
if Langevin == 1 
   gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
   derF = gradloglikHandle(model.Likelihood, model.y, F);
  
   % We assume that we precondition with K=Q^(-1)
   if strcmp(model.precondition, 'yes') == 1 
       der = derF'*model.K - F;
   else
       der = derF' - QF;
   end
   
   mu = F + (model.delta/2)*der;
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

for it = 1:(BurnInIters + Iters) 
%
    % STEP 1: sample
    if Langevin == 1
        % MALA  
        if strcmp(model.precondition, 'yes') == 1 
           Fnew = F + (model.delta/2)*der + sqrt(model.delta)*(model.Khalf*randn(n,1))';
        else
           Fnew = F + (model.delta/2)*der + sqrt(model.delta)*randn(1,n);
        end
    else
        % simple Metropolis proposal 
        Fnew = F + sqrt(model.delta).*randn(1,n);
    end
     
    % STEP 2: M-H step
   
    % perform an evaluation of the likelihood p(Y | F) 
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
    newLogLik = sum(newLogLik(:));  
    QFnew = Fnew*model.Q;
    QF2new = QFnew*Fnew';
    newLog_joint = newLogLik - 0.5*QF2new;

    % Metropolis-Hastings to accept-reject the proposal
    corrFactor = 0;
    if Langevin == 1 
        
        % new gradient 
        derFnew = gradloglikHandle(model.Likelihood, model.y, Fnew);
        
        % We assume that we precondition with K=Q^(-1)
        if strcmp(model.precondition, 'yes') == 1 
           dernew = derFnew'*model.K - Fnew;
        else
           dernew = derFnew' - QFnew; 
        end
        
        munew = Fnew + (model.delta/2)*dernew;
        
        %% compute the langevin proposal terms in the M-H ratio
        tmp = Fnew - mu; 
        if strcmp(model.precondition, 'yes') == 1
            tmp = tmp*model.Qhalf;
        end
        proposal_new_given_old = - (0.5/model.delta)*(tmp*tmp');
        tmp = F - munew;
        if strcmp(model.precondition, 'yes') == 1
            tmp = tmp*model.Qhalf;
        end
        proposal_old_given_new = - (0.5/model.delta)*(tmp*tmp');
        
        corrFactor = proposal_old_given_new - proposal_new_given_old; 
    end
    [accept, uprob] = metropolisHastings(newLog_joint+corrFactor, oldLog_joint, 0, 0); 
   
        
    acceptHistF(it) = accept;   
      
    if accept == 1
         F = Fnew;
         der = dernew; 
         mu = munew; 
         oldLogLik = newLogLik; 
         oldLog_joint = newLog_joint;
    end
    
    % Adapt proposal during burnin
    if mod(it,5) == 0 
        if (it >= 50) 
        accRateF = mean(acceptHistF((it-49):it))*100;
        if (it <= BurnInIters)
        if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
        %    
            model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;
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

