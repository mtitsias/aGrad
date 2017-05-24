
% MAIN FUNCTIONS THAT IMPLEMENT THE SAMPLING ALGORITHMS 
gpsampAuxZ_fixedhypers.m  % it implements aGrad-z for fixed kernel hyperparameters 
gpsampAuxU_fixedhypers.m  % it implements aGrad-u for fixed kernel hyperparameters 
gpsampAuxMarg_fixedhypers.m  % it implements mGrad for fixed kernel hyperparameters 
gpsampMALA_fixedhypers.m     % it implements pMALA for fixed kernel hyperparameters 
gpsampEllipt_fixedhypers.m   % it implements Ellipt for fixed kernel hyperparameters 
gpsamppCN_fixedhypers.m      % it implements both pCN (with the option langevin=0; see demos) 
                             % and pCNL (with the option langevin=1; see demos) for fixed kernel hyperparameters 
gpsampAuxZ_learnhypers_joint.m  % it implements aGrad-z-joint for learning hyperparameters
gpsampAuxZ_learnhypers.m        % it implements aGrad-z-gibbs for learning hyperparameters
gpsamppCN_learnhypers.m         % it implements pCN(L)-gibbs for learning hyperparameters
gpsampEllipt_learnhypers.m      % it implements Ellipt-gibbs for learning hyperparameters
gpsampAuxU_learnhypers.m        % it implements uGrad-z-gibbs for learning hyperparameters (this method was not used in the paper)

% DEMOS 
demos_informlik.m  % Runs all demos from section 5.1 Gaussian process regression and small noise limit
    % the individual demos are the following
    demRegressInformLikelMarg_fixedhypers.m
    demRegressInformLikelAuxZ_fixedhypers.m
    demRegressInformLikelAuxU_fixedhypers.m
    demRegressInformLikelMALA_fixedhypers.m
    demRegressInformLikelEllipt_fixedhypers.m
    demRegressInformLikelpCN_fixedhypers.m
    demRegressInformLikelpCNL_fixedhypers.m
demos_girolami.m  % Runs all demos from section 5.2 Log-Gaussian Cox Process (apart from the mMALA and RMHMC)
                  % for which you need to use Girolami and Calderhead software) 
    % the individual demos are the following            
    demLogGaussianCoxGirolamiMarg.m
    demLogGaussianCoxGirolamiAuxU.m
    demLogGaussianCoxGirolamiAuxZ.m
    demLogGaussianCoxGirolamiMALA.m
    demLogGaussianCoxGirolamiEllipt.m
    demLogGaussianCoxGirolamipCN.m
    demLogGaussianCoxGirolamipCNL.m 
demos_binaryclassification_fixedhypers.m % Runs all demos from section 5.3 Binary Gaussian process classification
    % the individual demos are the following
    demBinaryClassification_MALA_fixedhypers.m
    demBinaryClassification_Marg_fixedhypers.m
    demBinaryClassification_AuxZ_fixedhypers.m
    demBinaryClassification_AuxU_fixedhypers.m
    demBinaryClassification_Ellipt_fixedhypers.m
    demBinaryClassification_pCN_fixedhypers.m
    demBinaryClassification_pCNL_fixedhypers.m
demos_mnistsoftmax.m % Runs all MNIST demos from section 5.4 Sampling hyperparameters for binary 
                     % and multi-class Gaussian process classification
    demMultiClassMnist_AuxZadvanced_learnhypers.m
    demMultiClassMnist_Ellipt_learnhypers.m
    demMultiClassMnist_AuxZ_learnhypers.m 
    demMultiClassMnist_pCNL_learnhypers.m
demos_binaryclassification_learnhypers.m % runs demos for learing hyperparameter the five binary classfication datasets
                                         % only the for Pima datastes are discussed in the paper (section 5.4) 
    demBinaryClassification_AuxZ_learnhypers.m
    demBinaryClassification_AuxU_learnhypers.m
    demBinaryClassification_Ellipt_learnhypers.m
    demBinaryClassification_pCN_learnhypers.m
    demBinaryClassification_pCNL_learnhypers.m
    demBinaryClassification_AuxZadvanced_learnhypers.m
make_plots_informlik.m % create the plots and tables in section 5.1 (stored inside diagrams)     
make_plots_girolami.m  % creates the plots and tables in section 5.2 
make_plots_binaryclassification_fixedhypers.m  % creates the plots and tables in section 5.3 
make_plots_mnistsoftmax_learnhypers.m % creates the plot for mnist dataset in section 5.4 
make_plots_binaryclassification_learnhypers.m  % creates the plots and tables for the case we learn hyperparameters 
                                               % for all binary classification datasets 
                                               % (only the ones for the dataset Pima are mentioned in the paper)
make_acf.m % it creates the autocorrelation plots against CPU time that appear in Figure 1 in the main paper 

