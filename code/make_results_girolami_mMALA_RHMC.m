clear; 
outdir = '../diagrams/';
addpath ../results/; 
addpath ../diagrams/; 
addpath toolbox;
 
datasetName = 'GaussianCox';

Repeats = 10;

for i=1:2
%
  for rep=1:Repeats
   
  if i == 2    
     load(['RMHMC_LV_LogCox_reduced_repeat' num2str(rep) '.mat']);
  elseif i == 1
     load(['RMHMC_LV_LogCox_repeat' num2str(rep) '.mat']);
  end
  samples.F = xSaved; 

  samples.LogL =  histoldL;

  mcmcoptions.T = 5000;
  mcmcoptions.Burnin = 2000;
  mcmcoptions.StoreEvery = 1;
  
  if ~exist(['../results/logGaussianCoxGirolami_RHMCAndmMALA_' num2str(i) '.mat'])
  % compute statistics 
  summaryRHMC{rep} = summaryStatistics(samples);
  summaryRHMC{rep}.elapsed = TimeTaken;
  summaryRHMC{rep}.delta = StepSize; 
  summaryRHMC{rep}.eff_LogL = mcmc_ess(histoldL(mcmcoptions.Burnin+1:end));
  else
  load(['../results/logGaussianCoxGirolami_RHMCAndmMALA_' num2str(i) '.mat']);
  end
  
  ESSminRHMC(rep) = min(summaryRHMC{rep}.eff_F); 
  ESSmedianRHMC(rep) = median(summaryRHMC{rep}.eff_F); 
  ESSmaxRHMC(rep) = max(summaryRHMC{rep}.eff_F); 
  ESSlogLRHMC(rep) = summaryRHMC{rep}.eff_LogL; 
  TrainTimeRHMC(rep) = summaryRHMC{rep}.elapsed;
  deltaRHMC(rep) = summaryRHMC{rep}.delta; 
  
  if i == 2    
     load(['mMALA_LV_LogCox_reduced_repeat' num2str(rep) '.mat']);
  elseif i == 1
     load(['mMALA_LV_LogCox_repeat' num2str(rep) '.mat']);
  end

  samples.F = xSaved; 
  samples.LogL =  histoldL;

  if ~exist(['../results/logGaussianCoxGirolami_RHMCAndmMALA_' num2str(i) '.mat'])
  % compute statistics 
  summarymMALA{rep} = summaryStatistics(samples);
  summarymMALA{rep}.elapsed = TimeTaken;
  summarymMALA{rep}.delta = StepSize; 
  summarymMALA{rep}.eff_LogL = mcmc_ess(histoldL(mcmcoptions.Burnin+1:end));
  end
 
  ESSminmMALA(rep) = min(summarymMALA{rep}.eff_F); 
  ESSmedianmMALA(rep) = median(summarymMALA{rep}.eff_F); 
  ESSmaxmMALA(rep) = max(summarymMALA{rep}.eff_F); 
  ESSlogmMALA(rep) = summarymMALA{rep}.eff_LogL; 
  TrainTimemMALA(rep) = summarymMALA{rep}.elapsed;
  deltamMALA(rep) = summarymMALA{rep}.delta; 
  
  end
  
  if ~exist(['../results/logGaussianCoxGirolami_RHMCAndmMALA_' num2str(i) '.mat'])
   save(['../results/logGaussianCoxGirolami_RHMCAndmMALA_' num2str(i) '.mat'], 'summarymMALA', 'summaryRHMC');
  end

  fid = fopen([outdir datasetName '_RHMCAndmMALA_table' num2str(i) '.txt'],'w'); 
  fprintf(fid,'Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.) \\\\ \n'); 
  fprintf(fid,'\\midrule\n');
  fprintf(fid,'mMALA  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f) \\\\ \n', ...
               mean(TrainTimemMALA), mean(deltamMALA), mean(ESSminmMALA), mean(ESSmedianmMALA), mean(ESSmaxmMALA), mean(ESSminmMALA./TrainTimemMALA), std(ESSminmMALA./TrainTimemMALA));
  %fprintf(fid,'   &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
  %             std(TrainTimemMALA), std(deltamMALA), std(ESSminmMALA), std(ESSmedianmMALA), std(ESSmaxmMALA), std(ESSminmMALA./TrainTimemMALA));            
  fprintf(fid,'RHMC  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f) \\\\ \n', ...
               mean(TrainTimeRHMC), mean(deltaRHMC), mean(ESSminRHMC), mean(ESSmedianRHMC), mean(ESSmaxRHMC), mean(ESSminRHMC./TrainTimeRHMC), std(ESSminRHMC./TrainTimeRHMC));   
  %fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
  %             std(TrainTimeRHMC), std(deltaRHMC), std(ESSminRHMC), std(ESSmedianRHMC), std(ESSmaxRHMC), std(ESSminRHMC./TrainTimeRHMC));         
  fclose(fid);  
%  
end
