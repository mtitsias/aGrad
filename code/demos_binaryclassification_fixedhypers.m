% - Binary Classification examples
%
randn('seed', 1);
rand('seed', 1);
% how many times to repeat the experiments 
Repeats = 10;
for r=1:Repeats
   demBinaryClassification_MALA_fixedhypers(r);
   demBinaryClassification_Marg_fixedhypers(r);
   demBinaryClassification_AuxZ_fixedhypers(r);
   demBinaryClassification_AuxU_fixedhypers(r);
   demBinaryClassification_Ellipt_fixedhypers(r);
   demBinaryClassification_pCN_fixedhypers(r); 
   demBinaryClassification_pCNL_fixedhypers(r); 
end
