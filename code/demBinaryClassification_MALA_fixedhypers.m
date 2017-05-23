function demBinaryClassification_MALA_fixedhypers(rep)

addpath toolbox/; 
addpath ../data/;
%addpath ../../../Software/RHMC_GirolamiCalderhead/Bayes_Log_Reg/MCMC/Data/
precondition = 'yes'; % 'yes' or 'no'

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
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernAustralian;
        
   case 'German'

        %Load and prepare train & test data
        load('german.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernGerman; 
         
   case 'Heart'
        
        %Load and prepare train & test data
        load('heart.mat');
        t=X(:,end);
        X(:,end)=[];

        % German Credit - replace all 1s in t with 0s
        t(find(t==1)) = 0;
        % German Credit - replace all 2s in t with 1s
        t(find(t==2)) = 1;  
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernHeart;
        
   case 'Pima'
       
        %Load and prepare train & test data
        load('pima.mat');
        t=X(:,end);
        X(:,end)=[];  
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernPima;
        
   case 'Ripley'

        %Load and prepare train & test data
        load('ripley.mat');
        t=X(:,end);
        X(:,end)=[];
        
        load('initkern_binaryclass.mat'); 
        initkern = initkernRipley;
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
   
   model = gpsampCreate(Y, X, options); 
   model.GP.jitter = 1e-8;
   %model.GP{j}.logtheta = [2*log((max(X) - min(X))*0.8) 0];
   model.GP.logtheta = initkern; %[2*mean(log((max(X) - min(X))*0.8)) 0];
   
   model.XX = -2*(X*X') + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   
   model.constraints.kernHyper = 'free';
   model.constraints.likHyper = 'free';
   
   mcmcoptions.T = 5000;
   mcmcoptions.Burnin = 5000;
   mcmcoptions.StoreEvery = 1;
   mcmcoptions.Langevin = 1;

   
   model.precondition = precondition; 
   % precompute the inverse covariacne matrix Q 
   model.K = kernCompute(model.GP, model.X);
   [model.U, model.Lambda, tmp] = svd(model.K);
   model.Lambda = diag(model.Lambda);
   model.Lambda(model.Lambda<1e-15) = 1e-15;
   if strcmp(model.precondition,'yes') == 1
      model.Khalf = model.U*diag(model.Lambda.^0.5); % useful for preconditioned MALA
      model.Qhalf = model.U*diag(1./(model.Lambda.^0.5)); % useful for preconditioned MALA
   end
   model.invLambda = 1./model.Lambda;
   model.Q = model.U*diag(model.invLambda)*model.U'; 

   tic;
   [model samples accRates] = gpsampMALA_fixedhypers(model, mcmcoptions);
   elapsedTime = toc;

   % compute statistics 
   summaryMALA = summaryStatistics(samples);
   summaryMALA.elapsed = elapsedTime;
   summaryMALA.accRates = accRates;
   summaryMALA.delta = model.delta; 
   summaryMALA.eff_LogL = mcmc_ess(samples.LogL(mcmcoptions.Burnin+1:end));

   save(['../results/' dataName{1} '_repeat' num2str(rep) '_MALA_fixedhypers.mat'], 'summaryMALA');
%    
end