outdir = '../diagrams/';
fontsz = 26;
addpath ../results/; 

Repeats = 10;

for datasetName = {'Australian' 'German' 'Heart' 'Pima' 'Ripley'}% 'Caravan'}
%

if strcmp(datasetName,'Ripley')
   ax = [1 10000 -90 -60]; 
elseif strcmp(datasetName,'Pima')
   ax = [1 10000 -260 -210]; 
elseif strcmp(datasetName,'Heart')
   ax = [1 10000 -120 -60]; 
elseif strcmp(datasetName,'German')
   ax = [1 10000 -500 -400];
elseif strcmp(datasetName,'Australian')
   ax = [1 10000 -250 -150];
end

for rep=1:Repeats
    
% fname = [datasetName{1} '_AuxSimple_fixedhypers.mat'];
% load(fname);
% ESSminAuxSimple = min(summaryAuxSimple.eff_F); 
% ESSmedianAuxSimple = median(summaryAuxSimple.eff_F); 
% ESSmaxAuxSimple = max(summaryAuxSimple.eff_F); 
% ESSlogLAuxSimple = summaryAuxSimple.eff_LogL; 
% TrainTimeAuxSimple = summaryAuxSimple.elapsed;

fname = [datasetName{1} '_repeat' num2str(rep) '_AuxZ_fixedhypers.mat']; 
load(fname);
ESSminAuxZ(rep) = min(summaryAuxZ.eff_F); 
ESSmedianAuxZ(rep) = median(summaryAuxZ.eff_F); 
ESSmaxAuxZ(rep) = max(summaryAuxZ.eff_F); 
ESSlogLAuxZ(rep) = summaryAuxZ.eff_LogL; 
TrainTimeAuxZ(rep) = summaryAuxZ.elapsed;
deltaAuxZ(rep) = summaryAuxZ.delta; 

fname = [datasetName{1} '_repeat' num2str(rep) '_AuxU_fixedhypers.mat'];
load(fname);
ESSminAuxU(rep) = min(summaryAuxU.eff_F); 
ESSmedianAuxU(rep) = median(summaryAuxU.eff_F); 
ESSmaxAuxU(rep) = max(summaryAuxU.eff_F); 
ESSlogLAuxU(rep) = summaryAuxU.eff_LogL; 
TrainTimeAuxU(rep) = summaryAuxU.elapsed;
deltaAuxU(rep) = summaryAuxU.delta; 

fname = [datasetName{1} '_repeat' num2str(rep) '_Marg_fixedhypers.mat'];
load(fname);
ESSminMarg(rep) = min(summaryMarg.eff_F); 
ESSmedianMarg(rep) = median(summaryMarg.eff_F); 
ESSmaxMarg(rep) = max(summaryMarg.eff_F); 
ESSlogLMarg(rep) = summaryMarg.eff_LogL; 
TrainTimeMarg(rep) = summaryMarg.elapsed;
deltaMarg(rep) = summaryMarg.delta;

fname = [datasetName{1} '_repeat' num2str(rep) '_MALA_fixedhypers.mat'];
load(fname);
ESSminMALA(rep) = min(summaryMALA.eff_F); 
ESSmedianMALA(rep) = median(summaryMALA.eff_F); 
ESSmaxMALA(rep) = max(summaryMALA.eff_F); 
ESSlogLMALA(rep) = summaryMALA.eff_LogL; 
TrainTimeMALA(rep) = summaryMALA.elapsed;
deltaMALA(rep) = summaryMALA.delta;

fname = [datasetName{1} '_repeat' num2str(rep) '_Ellipt_fixedhypers.mat'];
load(fname);
ESSminEllipt(rep) = min(summaryEllipt.eff_F); 
ESSmedianEllipt(rep) = median(summaryEllipt.eff_F); 
ESSmaxEllipt(rep) = max(summaryEllipt.eff_F); 
ESSlogLEllipt(rep) = summaryEllipt.eff_LogL; 
TrainTimeEllipt(rep) = summaryEllipt.elapsed;

fname = [datasetName{1} '_repeat' num2str(rep) '_pCN_fixedhypers.mat'];
load(fname);
ESSminpCN(rep) = min(summarypCN.eff_F); 
ESSmedianpCN(rep) = median(summarypCN.eff_F); 
ESSmaxpCN(rep) = max(summarypCN.eff_F); 
ESSlogLpCN(rep) = summarypCN.eff_LogL; 
TrainTimepCN(rep) = summarypCN.elapsed;
deltapCN(rep) = summarypCN.delta;

fname = [datasetName{1} '_repeat' num2str(rep) '_pCNL_fixedhypers.mat'];
load(fname);
ESSminpCNL(rep) = min(summarypCNL.eff_F); 
ESSmedianpCNL(rep) = median(summarypCNL.eff_F); 
ESSmaxpCNL(rep) = max(summarypCNL.eff_F); 
ESSlogLpCNL(rep) = summarypCNL.eff_LogL; 
TrainTimepCNL(rep) = summarypCNL.elapsed;
deltapCNL(rep) = summarypCNL.delta;



if rep == 1
    
figure; 
plot(summaryAuxZ.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('aGrad-z');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_auxZ_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_auxZ_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summaryAuxU.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('aGrad-u');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_auxU_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_auxU_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summaryMarg.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('mGrad');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_marg_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_marg_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summaryEllipt.LogL,'k'); 
xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('Ellipt');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_ellipt_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_ellipt_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summarypCN.LogL,'k'); 
xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('pCN');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_pCN_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_pCN_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summarypCNL.LogL,'k'); 
xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('pCNL');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_pCNL_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_pCNL_logL_fixedhypers.eps']);
system(cmd);

figure; 
plot(summaryMALA.LogL,'k'); 
xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('pMALA');
set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName{1} '_mala_logL_fixedhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_mala_logL_fixedhypers.eps']);
system(cmd);

end

end

%disp(datasetName{1});
%[min(ESSminAuxZ), max(ESSminAuxZ), mean(ESSminAuxZ), median(ESSminAuxZ)]
%[min(ESSminAuxU), max(ESSminAuxU), mean(ESSminAuxU), median(ESSminAuxU)]
%[min(ESSminMarg), max(ESSminMarg), mean(ESSminMarg), median(ESSminMarg)]

fid = fopen([outdir datasetName{1} '_table_fixedhypers.txt'],'w'); 
fprintf(fid,' Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.) \\\\ \n'); 
fprintf(fid,'\\midrule\n');
fprintf(fid,'aGrad-z  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeAuxZ), mean(deltaAuxZ), mean(ESSminAuxZ), mean(ESSmedianAuxZ), mean(ESSmaxAuxZ), mean(ESSminAuxZ./TrainTimeAuxZ), std(ESSminAuxZ./TrainTimeAuxZ));
%fprintf(fid,'   &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimeAuxZ), std(deltaAuxZ), std(ESSminAuxZ), std(ESSmedianAuxZ), std(ESSmaxAuxZ), std(ESSminAuxZ./TrainTimeAuxZ));            
fprintf(fid,'aGrad-u  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f) \\\\ \n', ...
            mean(TrainTimeAuxU), mean(deltaAuxU), mean(ESSminAuxU), mean(ESSmedianAuxU), mean(ESSmaxAuxU), mean(ESSminAuxU./TrainTimeAuxU), std(ESSminAuxU./TrainTimeAuxU));   
%fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%           std(TrainTimeAuxU), std(deltaAuxU), std(ESSminAuxU), std(ESSmedianAuxU), std(ESSmaxAuxU), std(ESSminAuxU./TrainTimeAuxU));         
fprintf(fid,'mGrad  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeMarg),  mean(deltaMarg), mean(ESSminMarg), mean(ESSmedianMarg), mean(ESSmaxMarg), mean(ESSminMarg./TrainTimeMarg), std(ESSminMarg./TrainTimeMarg));
%fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimeMarg),  std(deltaMarg), std(ESSminMarg), std(ESSmedianMarg), std(ESSmaxMarg), std(ESSminMarg./TrainTimeMarg));
fprintf(fid,'pMALA  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeMALA),  mean(deltaMALA), mean(ESSminMALA), mean(ESSmedianMALA), mean(ESSmaxMALA), mean(ESSminMALA./TrainTimeMALA), std(ESSminMALA./TrainTimeMALA));  
%fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimeMALA),  std(deltaMALA), std(ESSminMALA), std(ESSmedianMALA), std(ESSmaxMALA), std(ESSminMALA./TrainTimeMALA));  
fprintf(fid,'Ellipt  &  %1.1f  &   &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimeEllipt),   mean(ESSminEllipt), mean(ESSmedianEllipt), mean(ESSmaxEllipt), mean(ESSminEllipt./TrainTimeEllipt), std(ESSminEllipt./TrainTimeEllipt)); 
%fprintf(fid,'  &  %1.1f  &   &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimeEllipt),  std(ESSminEllipt), std(ESSmedianEllipt), std(ESSmaxEllipt), std(ESSminEllipt./TrainTimeEllipt)); 
fprintf(fid,'pCN  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimepCN), mean(deltapCN), mean(ESSminpCN), mean(ESSmedianpCN), mean(ESSmaxpCN), mean(ESSminpCN./TrainTimepCN), std(ESSminpCN./TrainTimepCN)); 
%fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimepCN), std(deltapCN), std(ESSminpCN), std(ESSmedianpCN), std(ESSmaxpCN), std(ESSminpCN./TrainTimepCN));     
fprintf(fid,'pCNL  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f (%1.2f)\\\\ \n', ...
            mean(TrainTimepCNL),  mean(deltapCNL), mean(ESSminpCNL), mean(ESSmedianpCNL), mean(ESSmaxpCNL), mean(ESSminpCNL./TrainTimepCNL), std(ESSminpCNL./TrainTimepCNL));
%fprintf(fid,'  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            std(TrainTimepCNL),  std(deltapCNL), std(ESSminpCNL), std(ESSmedianpCNL), std(ESSmaxpCNL), std(ESSminpCNL./TrainTimepCNL));
fclose(fid);  


%
end

