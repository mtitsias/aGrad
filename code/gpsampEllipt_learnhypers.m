function [model samples accRates] = gpsampEllipt_learnhypers(model, trainOps)
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

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F');
oldLogLik = sum(oldLogLik(:)); 


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
   for j=1:model.J
     oldLogPriorLik(j) = sum( feval(lnpriorLik, model.Likelihood.logtheta, model.prior.likParams.a, model.prior.likParams.b) );      
   end
end     


cnt = 0;


[n D] = size(model.X);


for j=1:model.J
   GPtmp = model.GP{j};
   if isfield(model,'XX')
     GPtmp.XX = model.XX;
   end
   model.GP{j}.K = kernCompute(GPtmp, model.X);

   [L,er]=chol(model.GP{j}.K);
   model.GP{j}.L = L;
   L = L';
   % evaluate the new log GP prior value 
   model.GP{j}.invL = L\eye(n);
   model.GP{j}.LogDetK = 2*sum(log(diag(L)));
end 


% proposal for the kernel hyperparameters
model.kerndelta = 0.2*(1/model.prior.kernParams.b)*ones(1,model.J);
if model.Likelihood.nParams > 0
  model.likdelta = 0.2*(1/model.prior.likParams.b);
end

acceptHistLik = zeros(1, BurnInIters + Iters);
acceptHistKern = zeros(model.J, BurnInIters + Iters);
accRateK = ones(1,model.J);

range = 0.05; 

optLik = 0.25; 
optKern = 0.25;

nu = zeros(model.J,model.n);
angle_range = 0;
% adaption step size
epsilon = 0.05;

for it = 1:(BurnInIters + Iters) 
%
    % the code here follows Iain Murray's implementation 
    for j=1:model.J 
        nu(j,:) = randn(1, model.n)*model.GP{j}.L;
    end
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
        oldLogLik = loglikHandle(model.Likelihood, Y, Fnew');
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
        %
       end
       acceptHistLik(it) = accept; 
       %
    end
    
    % sample kernel hyperparameters
    if strcmp(model.constraints.kernHyper, 'free')
       %     
       GPnew = model.GP;  
       oldlogGP = zeros(1,model.J);
       newlogGP = zeros(1,model.J); 
       newLogPriorK = zeros(1,model.J); 
       for j=1:model.J
          GPnew{j}.logtheta = model.GP{j}.logtheta  + sqrt(model.kerndelta(j)).*randn(1,model.GP{j}.nParams);
          GPtmp = GPnew{j};
          if isfield(model,'XX')
            GPtmp.XX = model.XX;
          end
          GPnew{j}.K = kernCompute(GPtmp, X);   
         
          newL=jitterChol(GPnew{j}.K);
          GPnew{j}.L = newL;
          newL = newL';
          % evaluate the new log GP prior value 
          GPnew{j}.invL = newL\eye(n);
          GPnew{j}.LogDetK = 2*sum(log(diag(newL)));
      
          newlogGP(j) = - 0.5*GPnew{j}.LogDetK;
          oldlogGP(j) = - 0.5*model.GP{j}.LogDetK;
          temp = GPnew{j}.invL*F(j,:)'; 
          newlogGP(j) = newlogGP(j) - 0.5*(temp'*temp);
          temp = model.GP{j}.invL*F(j,:)'; 
          oldlogGP(j) = oldlogGP(j) - 0.5*(temp'*temp); 
          newLogPriorK(j) = sum( feval(lnpriorK, GPnew{j}.logtheta, model.prior.kernParams.a, model.prior.kernParams.b) );
          
          % Metropolis-Hastings to accept-reject the proposal
          accept = metropolisHastings(newlogGP(j) + newLogPriorK(j), oldlogGP(j) + oldLogPriorK(j), 0, 0);
          if accept == 1
          % 
             model.GP{j} = GPnew{j};
             oldLogPriorK(j) = newLogPriorK(j);  
          %     
          end
          acceptHistKern(j,it) = accept; 
       %
       end
    %
    end
  
   
    if mod(it,5) == 0
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
    
    % keep samples after burn in
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
          fprintf(1,'Iters=%d, AccRate for kernel hypers=%f, logL=%f\n',it, min(accRateK), oldLogLik);
       else
          printf(1,'Iters=%d, logL=%f\n',it, oldLogLik);
       end
    end
%    
end
%

model.F = F;

if strcmp(model.constraints.kernHyper, 'free')
   accRates.kern = mean(acceptHistKern(:,BurnInIters+1:end),2)*100;
   accRates.kern = accRates.kern';
end
if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
   accRates.lik = mean(acceptHistLik(BurnInIters+1:end))*100; 
end

model.acceptHistKern = acceptHistKern; 
model.acceptHistLik = acceptHistLik;

