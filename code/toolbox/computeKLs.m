function KLdiv = computeKLs(model, samples)
% auxiliary function
%
%
%

X = model.X;
Y = model.y;
  
[N D] = size(X);
  
sigma2n = exp(model.Likelihood.logtheta(1));
sigmaf = exp(model.GP.logtheta(D+1));

% GP posterior 
Knn = kernCompute(model.GP, X, X);

jitter = 0.000001*mean(diag(Knn));

L = chol(Knn + sigma2n*eye(N))';
alpha = L'\(L\Y);
muGP = Knn*alpha;   
v = L\Knn;  
CovarGP = Knn - v'*v + sigma2n*eye(N); 
     
Fs = samples.F;
muMCMC = mean(Fs);
CovarMCMC = cov(Fs,1)+sigma2n*eye(N);
    
KLdiv = kl(muGP, CovarGP, muMCMC, CovarMCMC);

