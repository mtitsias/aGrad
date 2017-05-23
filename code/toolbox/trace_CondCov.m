function [f, df] = trace_CondCov(Xu, X, logtheta)
% Description:  Function and gradient evalaution with respect to  pseudo-inputs 
%               of the trace of the conditional prior (given inducing variables)
%               for the expenential kernel with varied length-scale (ARD kernel)
%

% number of examples and dimension of input space
[n, D] = size(X);
m = round(size(Xu,1)/D);
jitter = 1e-8;

Xu = reshape(Xu,m,D);

sigmaf = exp(2*logtheta(D+1));
X = X ./ repmat(exp(logtheta(1:D)),n,1);
Xu = Xu ./ repmat(exp(logtheta(1:D)),m,1);

Kmm = Xu*Xu';
Kmm = repmat(diag(Kmm),1,m) + repmat(diag(Kmm)',m,1) - 2*Kmm;
Kmm = sigmaf*exp(-0.5*Kmm);
Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
Knm = sigmaf*exp(-0.5*Knm');


[Lkmm er] = jitterChol(Kmm);
%if er > 0 % add jitter
%   Kmm = Kmm + 1e-7*sigmaf*eye(m);
%   disp('jitter added')
%   Lkmm = chol(Kmm);
%end
%

Cnm1 = Knm/Lkmm;
Cmnmn = Cnm1'*Cnm1;

% value of the objective function
f = n*sigmaf - sum(diag(Cmnmn));

% compute derivatives
Pmnmn = (Lkmm\Cmnmn)';
BB1 = Lkmm\Pmnmn; 
BB1 = Kmm.*BB1;

Cnm1 = (Lkmm\Cnm1')';
Cnm1 = Cnm1.*Knm;
%
for d=1:D
    %
    % pseudo inputs derivatives
    Knm = -((repmat(Xu(:,d)',n,1)-repmat(X(:,d),1,m))/exp(logtheta(d)));
    Kmm = -((repmat(Xu(:,d)',m,1)-repmat(Xu(:,d),1,m))/exp(logtheta(d)));       
    
    dXu(:,d) = -(sum(Knm.*Cnm1,1) - sum(Kmm.*BB1,1))'; 
    %
end
%
dXu = 2*dXu;
df = reshape(dXu, m*D, 1);
