% -- Regression with different levels of noise 
%
randn('seed', 1);
rand('seed', 1);
% how many times to repeat the experiments 
Repeats = 10;
for r=1:Repeats
   demRegressInformLikelMarg_fixedhypers(r);
   demRegressInformLikelAuxZ_fixedhypers(r);
   demRegressInformLikelAuxU_fixedhypers(r);
   demRegressInformLikelMALA_fixedhypers(r);
   demRegressInformLikelEllipt_fixedhypers(r);
   demRegressInformLikelpCN_fixedhypers(r);
   demRegressInformLikelpCNL_fixedhypers(r);
end
