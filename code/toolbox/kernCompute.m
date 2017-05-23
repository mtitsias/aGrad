function Knm = covfuncCompute(GPprior, X, Xu)
%function Knm = covfuncCompute(logtheta, X, Xu)
%
%Description:  It computes the covariance function between 
%              two set of inputs points: X and Xu.  
%
%Supported covariance functions:  RBF and ARD kernel.
%      hyperparameters are modelled in the log space  
%      hyp(1:D) = log( ell^2 )
%      hyp(D+1) = log( sigmaf2 )
%
%      where cov = sigmaf2 * exp[-0.5 * sum_i (x_i - x'_i)^2/ell^2]

 

jitter = GPprior.jitter; 

switch GPprior.type 
  case 'seard'
    
    [n D] = size(X);
    logtheta = GPprior.logtheta(:);
    sigmaf = exp(logtheta(D+1));
    X = X ./ repmat(exp(logtheta(1:D)/2)',n,1);
    
    if nargin == 3
      [m,D] = size(Xu);   
      Xu = Xu ./ repmat(exp(logtheta(1:D)/2)',m,1);
      %
      Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
      Knm = sigmaf*exp(-0.5*Knm');
    else
      Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
      Knm = sigmaf*exp(-0.5*Knm') + jitter*eye(n); 
    end
 case 'se'
      
    [n D] = size(X);
    logtheta = GPprior.logtheta(:);
    sigmaf = exp(logtheta(2));
    
    if nargin == 3
      [m,D] = size(Xu);   
      Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
      Knm = sigmaf*exp(-Knm'./(2*exp(logtheta(1))) );
    else
      if isfield(GPprior,'XX') 
          Knm = sigmaf*exp(-GPprior.XX./(2*exp(logtheta(1))) ) + jitter*eye(n);
      else
          Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
          Knm = sigmaf*exp(-Knm'./(2*exp(logtheta(1))) ) + jitter*eye(n);
      end
    end
  case 'sre'
      
    [n D] = size(X);
    logtheta = GPprior.logtheta(:);
    sigmaf = exp(logtheta(2));
    
    if nargin == 3
      [m,D] = size(Xu);   
      Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
      Knm = Knm.^0.5;
      Knm = sigmaf*exp(-Knm'./(2*exp(logtheta(1))) );
    else
      Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
      Knm = Knm.^0.5;
      Knm = sigmaf*exp(-Knm'./(2*exp(logtheta(1))) ); 
    end
        
end  