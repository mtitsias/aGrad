close all; 
addpath toolbox/; 
addpath ../data/;

for dataName = {'Australian' 'German' 'Heart' 'Pima' 'Ripley'}% 'Caravan'}
%     
   disp(['Running with ' dataName{1} ' dataset']);
   switch(dataName{1})
    
   case 'Caravan'
    
        %Load and prepare train & test data
        load('caravan.mat');
        
        % End column contains binary response data
        t=X(:,end);
        X(:,end) = [];

   case 'Australian'
    
        %Load and prepare train & test data
        load('australian.mat');
        t=X(:,end);
        X(:,end)=[];
        
   case 'German'

        %Load and prepare train & test data
        load('german.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;
         
   case 'Heart'
        
        %Load and prepare train & test data
        load('heart.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;
        
   case 'Pima'
       
        %Load and prepare train & test data
        load('pima.mat');
        t=X(:,end);
        X(:,end)=[];
        
   case 'Ripley'

        %Load and prepare train & test data
        load('ripley.mat');
        t=X(:,end);
        X(:,end)=[];

   end   
   
   Y = 2*t - 1;
   
   [n D] = size(X);
   
   % normalize the data so that the inputs have unit variance and zero mean  
   meanX = mean(X);
   sqrtvarX = sqrt(var(X)); 
   X = X - repmat(meanX, n, 1);
   X = X./repmat(sqrtvarX, n, 1);
   
   % model options
   options = gpsampOptions('classification'); 
   options.kern = 'se';
   
   model = gpsampCreate(Y, X, options, 1); 
   for j=1:model.J 
     model.GP{j}.jitter = 1e-8;
     %model.GP{j}.logtheta = [2*log((max(X) - min(X))*0.8) 0];
     model.GP{j}.logtheta = [2*mean(log((max(X) - min(X))*0.8)) 0];
   end
   model.XX = -2*(X*X') + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   
   model.constraints.kernHyper = 'free';
   model.constraints.likHyper = 'free';
   
   mcmcoptions.T = 10000;
   mcmcoptions.Burnin = 40000;
   mcmcoptions.StoreEvery = 1;
   mcmcoptions.Langevin = 1;

   tic;
   [model samples accRates] = gpsampAuxU_learnhypers(model, mcmcoptions);
   elapsedTime = toc;

   % compute statistics  
   samples.F = squeeze(samples.F)';
   summaryAuxU = summaryStatistics(samples);
   summaryAuxU.elapsed = elapsedTime;
   summaryAuxU.accRates = accRates;
   summaryAuxU.delta = model.delta; 
   summaryAuxU.eff_LogL = mcmc_ess(samples.LogL(mcmcoptions.Burnin+1:end));
   summaryAuxU.kerndelta = model.kerndelta;
   for j=1:model.GP{1}.nParams
      summaryAuxU.eff_kern(j) = mcmc_ess(samples.kernLogtheta{1}(:,j));
   end
   
   save(['../results/' dataName{1} '_AuxU_learnhypers.mat'], 'summaryAuxU', 'model', 'samples');
%    
end