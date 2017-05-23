function [model samples accRates] = gpsampAuxU_learnhypers(model, trainOps)
% It samples also hyperparameters (Gibbs-based)


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F = zeros(model.J, n, num_stored);
for j=1:model.J
  samples.kernLogtheta{j} = zeros(num_stored, model.GP{j}.nParams);
end
samples.likLogtheta = zeros(num_stored, model.Likelihood.nParams);
samples.LogL = zeros(1, num_stored);
samples.ff = zeros(1, num_stored);
X = model.X;
Y = model.y;
F = model.F; 


U = zeros(size(F));
twodeltaU = zeros(size(F));
Fnew = zeros(size(F));
derFAderF = zeros(1, model.J);
derFAderFnew = zeros(1, model.J);
 
% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F');
oldLogLik = sum(oldLogLik(:)); 

Langevin = trainOps.Langevin;

oldLogPriorK = zeros(1,model.J);
% evaluation of the log prior for the kernel hyperparameters
if strcmp(model.constraints.kernHyper, 'free') 
   lnpriorK = ['ln', model.prior.kernParams.type,'pdf'];
   for j=1:model.J
     oldLogPriorK(j) = sum( feval(lnpriorK, model.GP{j}.logtheta, model.prior.kernParams.a, model.prior.kernParams.b) );
   end
end

oldLogPriorLik = zeros(1,model.J);
% evaluation of the log prior for the likelihood hyperparameters
if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
   lnpriorLik = ['ln', model.prior.likParams.type,'pdf'];
   oldLogPriorLik = feval(lnpriorLik, model.Likelihood.logtheta, model.prior.likParams.a, model.prior.likParams.b);      
end     


cnt = 0;

[n D] = size(model.X);

model.delta = 2/sqrt(n);

for j=1:model.J
   GPtmp = model.GP{j};
   if isfield(model,'XX')
     GPtmp.XX = model.XX;
   end
   model.GP{j}.K = kernCompute(GPtmp, model.X);
   
   L = chol(model.GP{j}.K + (model.delta/2)*eye(n));
   tmp = L\eye(n);
   tmp = tmp*(model.delta/2);
   model.GP{j}.invK = (model.delta/2)*eye(n) - tmp*tmp';
   model.GP{j}.auxPostL = jitterChol(model.GP{j}.invK);
   
   if strcmp(model.constraints.kernHyper, 'free')
      % the prior appears in the acceptacne ratio 
      [L,er]=jitterChol(model.GP{j}.K);  
      L = L';
      model.GP{j}.LogDetK = 2*sum(log(diag(L)));
      model.GP{j}.invL = L\eye(n);
   end
end 

% If the Langevin algorithm is used then 
% gradient wrt latent variables are required 
if Langevin == 1 
   gradloglikHandle = str2func(['grad', 'logL' model.Likelihood.type]);
   derF = gradloglikHandle(model.Likelihood, model.y, F'); 
   for j=1:model.J    
       tmp = model.GP{j}.auxPostL*derF(:,j);
       derFAderF(j) = 0.5*(tmp'*tmp);
   end       
end

% proposal for the kernel hyperparameters
model.kerndelta = 0.2*(1/model.prior.kernParams.b)*ones(1,model.J);
if model.Likelihood.nParams > 0
  model.likdelta = 0.2*(1/model.prior.likParams.b);
end


acceptHistF = zeros(1, BurnInIters + Iters);
acceptHistLik = zeros(1, BurnInIters + Iters);
acceptHistKern = zeros(model.J, BurnInIters + Iters);
accRateK = ones(1,model.J);

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

optLik = 0.25; 
optKern = 0.25;

% adaption step size
epsilon = 0.05;
 
for it = 1:(BurnInIters + Iters) 
%
    % STEP 1: Sample auxiliary variables
    for j=1:model.J  
        U(j,:) = F(j,:) + sqrt(model.delta/2)*randn(1,n);        
        twodeltaU(j,:) = (2/model.delta)*U(j,:);        
    end
    
    
    % STEP 2: M-H step
    for j=1:model.J 
       Fnew(j,:) = (randn(1, n) + (twodeltaU(j,:) + derF(:,j)')*model.GP{j}.auxPostL')*model.GP{j}.auxPostL;
    end
 
    % perform an evaluation of the likelihood p(Y | F) 
    newLogLik = loglikHandle(model.Likelihood, Y, Fnew');
    newLogLik = sum(newLogLik(:));  
    
    % Metropolis-Hastings to accept-reject the proposal
    corrFactor = 0;
    if Langevin == 1 
        
        % new gradient 
        derFnew = gradloglikHandle(model.Likelihood, model.y, Fnew);
          
        corrFactor = 0;
        for j=1:model.J    
          tmp = model.GP{j}.auxPostL*derFnew(:,j);
          derFAderFnew(j) = 0.5*(tmp'*tmp);
          tmp = (twodeltaU(j,:)*model.GP{j}.auxPostL')*model.GP{j}.auxPostL;
        
          rhoxyu = (F(j,:) - tmp)*derFnew(:,j) - derFAderFnew(j);
          rhoyxu = (Fnew(j,:) - tmp)*derF(:,j) - derFAderF(j);
        
          corrFactor = corrFactor + rhoxyu - rhoyxu;  
        end
    end
    [accept, uprob] = metropolisHastings(newLogLik+corrFactor, oldLogLik, 0, 0); 
    
   
    if accept == 1
         F = Fnew;
         oldLogLik = newLogLik; 
         derF = derFnew;
         derFAderF = derFAderFnew;
    end
    acceptHistF(it) = accept;     
    
    
    % sample gp likelihood parameters 
    if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
       % 
       newlogLik = randn(1,model.Likelihood.nParams).*sqrt(model.likdelta) + model.Likelihood.logtheta;
       
       Lik1 = model.Likelihood;
       Lik1.logtheta = newlogLik;
       % perform an evaluation of the likelihood p(Y | F) 
       newLogLik = loglikHandle(Lik1, Y, F'); 
       newL = newLogLik;
       newLogLik = sum(newLogLik(:)); 
       
       newLogPriorLik = feval(lnpriorLik, newlogLik, model.prior.likParams.a, model.prior.likParams.b);
       
       % Metropolis-Hastings to accept-reject the proposal
       oldlogP = oldLogLik + sum(oldLogPriorLik(:));
       newlogP = newLogLik + sum(newLogPriorLik(:)); 
       %
       [accept, uprob] = metropolisHastings(newlogP, oldlogP, 0, 0);
       if accept == 1
        %
          model.Likelihood.logtheta = newlogLik;
          oldLogPriorLik = newLogPriorLik;
          oldLogLik = newLogLik;
          if Langevin == 1          
             derF = gradloglikHandle(model.Likelihood, model.y, F');
             for j=1:model.J    
                tmp = model.GP{j}.auxPostL*derF(:,j);
                derFAderF(j) = 0.5*(tmp'*tmp);
             end    
          end
       end
       acceptHistLik(it) = accept; 
       %
    end
    
    % sample kernel hyperparameters
    if strcmp(model.constraints.kernHyper, 'free')
       % 
       Z = F;
       GPnew = model.GP;  
       oldlogGP = zeros(1,model.J);
       newlogGP = zeros(1,model.J); 
       newLogPriorK = zeros(1,model.J); 
       for j=1:model.J         
          GPnew{j}.logtheta = model.GP{j}.logtheta  + sqrt(model.kerndelta(j)).*randn(1,model.GP{j}.nParams);
          GPtmp = GPnew{j};
          GPtmp.XX = model.XX;
          GPnew{j}.K = kernCompute(GPtmp, X);   
          newL=jitterChol(GPnew{j}.K);
         
          newL = newL';
          % evaluate the new log GP prior value 
          GPnew{j}.invL = newL\eye(n);
          GPnew{j}.LogDetK = 2*sum(log(diag(newL)));
      
          newlogGP(j) = - 0.5*GPnew{j}.LogDetK;
          oldlogGP(j) = - 0.5*model.GP{j}.LogDetK;
          temp = GPnew{j}.invL*Z(j,:)'; 
          newlogGP(j) = newlogGP(j) - 0.5*(temp'*temp);
          temp = model.GP{j}.invL*Z(j,:)'; 
          oldlogGP(j) = oldlogGP(j) - 0.5*(temp'*temp); 
          newLogPriorK(j) = sum( feval(lnpriorK, GPnew{j}.logtheta, model.prior.kernParams.a, model.prior.kernParams.b) );
          
          % Metropolis-Hastings to accept-reject the proposal
          accept = metropolisHastings(newlogGP(j) + newLogPriorK(j), oldlogGP(j) + oldLogPriorK(j), 0, 0);
      
          %
          if accept == 1
          % 
             model.GP{j} = GPnew{j};
             oldLogPriorK(j) = newLogPriorK(j);    
             L = chol(model.GP{j}.K + (model.delta/2)*eye(n));
             tmp = L\eye(n);
             tmp = tmp*(model.delta/2);
             model.GP{j}.invK = (model.delta/2)*eye(n) - tmp*tmp';
             model.GP{j}.auxPostL = jitterChol(model.GP{j}.invK);
          %     
          end   
          acceptHistKern(j,it) = accept; 
       %
       end
    %
    end
  
    if mod(it,5) == 0 
    % Adapt the GP function proposal during burnin 
    if (it >= 50) 
        accRateF = mean(acceptHistF((it-49):it))*100;
        if (it <= BurnInIters)
        if (accRateF > (100*(opt+range))) || (accRateF < (100*(opt-range)))
            model.delta = model.delta + (epsilon*((accRateF/100 - opt)/opt))*model.delta;   
            for j=1:model.J
               L = chol(model.GP{j}.K + (model.delta/2)*eye(model.n));  
               tmp = L\eye(n);
               tmp = tmp*(model.delta/2);
               model.GP{j}.invK = (model.delta/2)*eye(n) - tmp*tmp';
               model.GP{j}.auxPostL = jitterChol(model.GP{j}.invK); 
            end
        end
        end
    end
    
    % Adapt the likelihoohd hypers proposal during burnin 
    if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
    if (it >= 50)
       accRateL = mean(acceptHistLik((it-49):it))*100;
       if (it <= BurnInIters)
          if (accRateL > (100*(optLik+range))) | (accRateL < (100*(optLik-range)))
             model.likdelta = model.likdelta + (epsilon*((accRateL/100 - optLik)/optLik))*model.likdelta;   
          end
       end
    end
    end  
   
    
    % Adapt the kernel hypers proposal during burnin 
    if strcmp(model.constraints.kernHyper, 'free')
    if (it >= 50)
       for j=1:model.J 
       accRateK(j) = mean(acceptHistKern(j,(it-49):it))*100;
       if (it <= BurnInIters)
       if (accRateK(j) > (100*(optKern+range))) || (accRateK(j) < (100*(optKern-range)))
           model.kerndelta(j) = model.kerndelta(j) + (epsilon*((accRateK(j)/100 - optKern)/optKern))*model.kerndelta(j);   
       end
       end
       end
    end
    end
    
    end
    
    
    % keep samples after burnin
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
    %
        cnt = cnt + 1; 
        samples.F(:,:,cnt) = F;
        if strcmp(model.constraints.kernHyper, 'free')
           %put all kernel hyperparameter in one vector 
           for j=1:model.J
             samples.kernLogtheta{j}(cnt,:) = model.GP{j}.logtheta; 
           end
        end 
        if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
           samples.likLogtheta(cnt,:) = model.Likelihood.logtheta;  
        end
    %
    end
    %  
    
    samples.LogL(it) = oldLogLik;       
    samples.ff(it) = 0.5*sum(sum(F.*F));
           
    if mod(it,500) == 0 
       if strcmp(model.constraints.kernHyper, 'free') 
          fprintf(1,'Iters=%d, AccRate for F=%f, AccRate for kernel hypers=%f, logL=%f\n',it, accRateF, min(accRateK), oldLogLik);
       else
          fprintf(1,'Iters=%d, AccRate for F=%f, logL=%f\n',it,accRateF, oldLogLik);
       end
    end    
%    
end
%

model.F = F;
accRates.F = mean(acceptHistF(BurnInIters+1:end))*100; 

if strcmp(model.constraints.kernHyper, 'free')
   accRates.kern = mean(acceptHistKern(:,BurnInIters+1:end),2)*100;
   accRates.kern = accRates.kern';
end
if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
   accRates.lik = mean(acceptHistLik(BurnInIters+1:end))*100; 
end

model.acceptHistF = acceptHistF;
model.acceptHistKern = acceptHistKern; 
model.acceptHistLik = acceptHistLik;

