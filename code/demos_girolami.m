% - Girolami log Gaussian Cox demos
%
randn('seed', 1);
rand('seed', 1);
% how many times to repeat the experiments 
Repeats = 10;
for r=1:Repeats
   demLogGaussianCoxGirolamiMarg(r);
   demLogGaussianCoxGirolamiAuxU(r);
   demLogGaussianCoxGirolamiAuxZ(r);
   demLogGaussianCoxGirolamiMALA(r);
   demLogGaussianCoxGirolamiEllipt(r);
   demLogGaussianCoxGirolamipCN(r); 
   demLogGaussianCoxGirolamipCNL(r); 
end
