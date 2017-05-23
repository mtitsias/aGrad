outdir = '../diagrams/';
fontsz = 26;
addpath ../results/; 

algnames = {'aGrad-z-gibbs', 'aGrad-z-joint', 'Ellipt-gibbs', 'pCNL-gibbs'};

datasetName = 'mnistsoftmax';
%

fname = [datasetName '_AuxZ_learnhypers.mat'];
load(fname);
fname = [datasetName '_AuxZadvanced_learnhypers.mat'];
load(fname);
fname = [datasetName '_Ellipt_learnhypers.mat'];
load(fname);
fname = [datasetName '_pCNL_learnhypers.mat'];
load(fname);
figure;
hold on; 
plot(summaryAuxZ.LogL,'k')
plot(summaryAuxZadvanced.LogL,'r');
plot(summaryEllipt.LogL,'b');
plot(summarypCNL.LogL,'g');
ylabel('Log-likelihood','fontsize',fontsz-18);
xlabel('Sampling iteration','fontsize',fontsz-18);
%legend('aGrad-z-gibbs','aGrad-z-joint','Ellipt', 'pCNL','Location','SouthEast');
axis([0 100000 -2000 0]); 
box on;
pbaspect([5 1 1]);
set(gca,'fontsize',fontsz-18);
%title('aGrad-z');
%set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName '_logLs_learnhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName '_logLs_learnhypers.eps']);
system(cmd);


kernlength  = [];
kernvar  = [];

%fname = [datasetName{1} '_AuxSimple_learnhypers.mat'];
%load(fname);
%ESSminAuxSimple = min(summaryAuxSimple.eff_F); 
%ESSmedianAuxSimple = median(summaryAuxSimple.eff_F); 
%ESSmaxAuxSimple = max(summaryAuxSimple.eff_F); 
%ESSlogLAuxSimple = summaryAuxSimple.eff_LogL; 
%TrainTimeAuxSimple = summaryAuxSimple.elapsed;

fname = [datasetName '_AuxZ_learnhypers.mat'];
load(fname);
ESSminAuxZ = min(summaryAuxZ.eff_F); 
ESSmedianAuxZ = median(summaryAuxZ.eff_F); 
ESSmaxAuxZ = max(summaryAuxZ.eff_F); 
%ESSlogLAuxZ = summaryAuxZ.eff_LogL; 
TrainTimeAuxZ = summaryAuxZ.elapsed;

kernlength = [kernlength, samples.kernLogtheta{1}(:,1)];
kernvar = [kernvar, samples.kernLogtheta{1}(:,2)];

fname = [datasetName '_AuxZadvanced_learnhypers.mat'];
load(fname);
ESSminAuxZadvanced = min(summaryAuxZadvanced.eff_F); 
ESSmedianAuxZadvanced = median(summaryAuxZadvanced.eff_F); 
ESSmaxAuxZadvanced = max(summaryAuxZadvanced.eff_F); 
%ESSlogLAuxZadvanced = summaryAuxZadvanced.eff_LogL; 
TrainTimeAuxZadvanced = summaryAuxZadvanced.elapsed;

kernlength = [kernlength, samples.kernLogtheta{1}(:,1)];
kernvar = [kernvar, samples.kernLogtheta{1}(:,2)];

%fname = [datasetName{1} '_AuxU_learnhypers.mat'];
%load(fname);
%ESSminAuxU = min(summaryAuxU.eff_F); 
%ESSmedianAuxU = median(summaryAuxU.eff_F); 
%ESSmaxAuxU = max(summaryAuxU.eff_F); 
%ESSlogLAuxU = summaryAuxU.eff_LogL; 
%TrainTimeAuxU = summaryAuxU.elapsed;

%fname = [datasetName{1} '_Marg_learnhypers.mat'];
%load(fname);
%ESSminMarg = min(summaryMarg.eff_F); 
%ESSmedianMarg = median(summaryMarg.eff_F); 
%ESSmaxMarg = max(summaryMarg.eff_F); 
%ESSlogLMarg = summaryMarg.eff_LogL; 
%TrainTimeMarg = summaryMarg.elapsed;

%fname = [datasetName{1} '_MALA_learnhypers.mat'];
%load(fname);
%ESSminMALA = min(summaryMALA.eff_F); 
%ESSmedianMALA = median(summaryMALA.eff_F); 
%ESSmaxMALA = max(summaryMALA.eff_F); 
%ESSlogLMALA = summaryMALA.eff_LogL; 
%TrainTimeMALA = summaryMALA.elapsed;

%fname = [datasetName '_Ellipt_learnhypers_betterInit.mat'];
%load(fname);
%ESSminEllipt = min(summaryEllipt.eff_F); 
%ESSmedianEllipt = median(summaryEllipt.eff_F); 
%ESSmaxEllipt = max(summaryEllipt.eff_F); 
%ESSlogLEllipt = summaryEllipt.eff_LogL; 
%TrainTimeEllipt = summaryEllipt.elapsed;

%kernlength = [kernlength, samples.kernLogtheta{1}(:,1)];
%kernvar = [kernvar, samples.kernLogtheta{1}(:,2)];


%fname = [datasetName{1} '_pCN_learnhypers.mat'];
%load(fname);
%ESSminpCN = min(summarypCN.eff_F); 
%ESSmedianpCN = median(summarypCN.eff_F); 
%ESSmaxpCN = max(summarypCN.eff_F); 
%ESSlogLpCN = summarypCN.eff_LogL; 
%TrainTimepCN = summarypCN.elapsed;

%fname = [datasetName '_pCNL_learnhypers_betterInit.mat'];
%load(fname);
%ESSminpCNL = min(summarypCNL.eff_F); 
%ESSmedianpCNL = median(summarypCNL.eff_F); 
%ESSmaxpCNL = max(summarypCNL.eff_F); 
%ESSlogLpCNL = summarypCNL.eff_LogL; 
%TrainTimepCNL = summarypCNL.elapsed;
%
%kernlength = [kernlength, samples.kernLogtheta{1}(:,1)];
%kernvar = [kernvar, samples.kernLogtheta{1}(:,2)];

fid = fopen([outdir datasetName '_table_learnhypers.txt'],'w'); 
fprintf(fid,' Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s \\\\ \n'); 
fprintf(fid,'\\midrule\n');
%fprintf(fid,'AMALA-simple  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeAuxSimple,  summaryAuxSimple.delta, ESSminAuxSimple, ESSmedianAuxSimple, ESSmaxAuxSimple, ESSminAuxSimple/TrainTimeAuxSimple);
fprintf(fid,'aGrad-z-gibbs  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
            TrainTimeAuxZ,  summaryAuxZ.delta, ESSminAuxZ, ESSmedianAuxZ, ESSmaxAuxZ, ESSminAuxZ/TrainTimeAuxZ);
fprintf(fid,'aGrad-z-joint  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
            TrainTimeAuxZadvanced,  summaryAuxZadvanced.delta, ESSminAuxZadvanced, ESSmedianAuxZadvanced, ESSmaxAuxZadvanced, ESSminAuxZadvanced/TrainTimeAuxZadvanced);       
%fprintf(fid,'mGrad  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeMarg,  summaryMarg.delta, ESSminMarg, ESSmedianMarg, ESSmaxMarg, ESSminMarg/TrainTimeMarg);
%fprintf(fid,'pMALA  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeMALA,  summaryMALA.delta, ESSminMALA, ESSmedianMALA, ESSmaxMALA, ESSminMALA/TrainTimeMALA);  
%       fprintf(fid,'Ellipt-gibbs  &  %1.1f   &    &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeEllipt, ESSminEllipt, ESSmedianEllipt, ESSmaxEllipt, ESSminEllipt/TrainTimeEllipt); 
%      fprintf(fid,'pCN  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimepCN,  summarypCN.delta, ESSminpCN, ESSmedianpCN, ESSmaxpCN, ESSminpCN/TrainTimepCN); 
%       fprintf(fid,'pCNL-gibbs  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimepCNL,  summarypCNL.delta, ESSminpCNL, ESSmedianpCNL, ESSmaxpCNL, ESSminpCNL/TrainTimepCNL);
fclose(fid);  


fid = fopen([outdir datasetName '_table_learnhypers_kern.txt'],'w'); 
fprintf(fid,' Method &  Time(s) & Step size $\\kappa$  &  ESS (Min, Med, Max)  & Min ESS/s \\\\ \n'); 
fprintf(fid,'\\midrule\n');
%fprintf(fid,'AMALA-simple  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeAuxSimple,  summaryAuxSimple.delta, ESSminAuxSimple, ESSmedianAuxSimple, ESSmaxAuxSimple, ESSminAuxSimple/TrainTimeAuxSimple);
fprintf(fid,'aGrad-z-gibbs  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
            TrainTimeAuxZ,  mean(summaryAuxZ.kerndelta), min(summaryAuxZ.eff_kern), min(summaryAuxZ.eff_kern), median(summaryAuxZ.eff_kern), max(summaryAuxZ.eff_kern)/TrainTimeAuxZ);
fprintf(fid,'aGrad-z-joint  &   %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
            TrainTimeAuxZadvanced,  mean(summaryAuxZadvanced.kerndelta), min(summaryAuxZadvanced.eff_kern), median(summaryAuxZadvanced.eff_kern), max(summaryAuxZadvanced.eff_kern), min(summaryAuxZadvanced.eff_kern)/TrainTimeAuxZadvanced);       
%fprintf(fid,'mGrad  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeMarg,  summaryMarg.delta, ESSminMarg, ESSmedianMarg, ESSmaxMarg, ESSminMarg/TrainTimeMarg);
%fprintf(fid,'pMALA  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeMALA,  summaryMALA.delta, ESSminMALA, ESSmedianMALA, ESSmaxMALA, ESSminMALA/TrainTimeMALA);  
%       fprintf(fid,'Ellipt-gibbs  &  %1.1f   &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimeEllipt, summaryEllipt.kerndelta,  min(summaryEllipt.eff_kern), median(summaryEllipt.eff_kern), max(summaryEllipt.eff_kern), min(summaryEllipt.eff_kern)/TrainTimeEllipt); 
%      fprintf(fid,'pCN  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimepCN,  summarypCN.delta, ESSminpCN, ESSmedianpCN, ESSmaxpCN, ESSminpCN/TrainTimepCN); 
%       fprintf(fid,'pCNL-gibbs  &  %1.1f  &  %1.3f  &  (%1.1f, %1.1f, %1.1f)  &  %1.2f \\\\ \n', ...
%            TrainTimepCNL,  summarypCNL.kerndelta, min(summarypCNL.eff_kern), median(summarypCNL.eff_kern), max(summarypCNL.eff_kern), min(summarypCNL.eff_kern)/TrainTimepCNL);
fclose(fid);  

ax = [0 20000 -300 -100]; 

% figure;
% boxplot(kernlength, algnames(1:2));
% set(findobj(gca,'Type','text'),'FontSize',fontsz); 
% xlabel('Algorithm','fontsize',fontsz);
% ylabel('log-lengthscale','fontsize',fontsz);
% set(gca,'fontsize',fontsz);
% %title('AuxLang-Z');
% print('-depsc2', '-r300', [outdir datasetName '_kernlength']);
% cmd = sprintf('epstopdf %s', [outdir datasetName '_kernlength.eps']);
% system(cmd);
% 
% figure;
% boxplot(kernvar, algnames(1:2));
% set(findobj(gca,'Type','text'),'FontSize',fontsz); 
% xlabel('Algorithm','fontsize',fontsz);
% ylabel('log-variance','fontsize',fontsz);
% set(gca,'fontsize',fontsz);
% %title('AuxLang-Z');
% print('-depsc2', '-r300', [outdir datasetName '_kernvar']);
% cmd = sprintf('epstopdf %s', [outdir datasetName '_kernvar.eps']);
% system(cmd);



figure; 
plot(summaryAuxZadvanced.LogL(end-20000:end),'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('aGrad-z-joint');
%set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName '_auxZadvanced_logL_learnhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName '_auxZadvanced_logL_learnhypers.eps']);
system(cmd);

figure; 
plot(summaryAuxZ.LogL(end-20000:end),'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log likelihood','fontsize',fontsz);
axis(ax);
set(gca,'fontsize',fontsz);
title('aGrad-z-gibbs');
%set(gca, 'XTick', [3000 6000 9000]);
print('-depsc2', '-r300', [outdir datasetName '_auxZ_logL_learnhypers']);
cmd = sprintf('epstopdf %s', [outdir datasetName '_auxZ_logL_learnhypers.eps']);
system(cmd);

% figure; 
% plot(summaryAuxU.LogL,'k'); 
% %xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% axis(ax);
% set(gca,'fontsize',fontsz);
% title('aGrad-u');
% set(gca, 'XTick', [3000 6000 9000]);
% print('-depsc2', '-r300', [outdir datasetName{1} '_auxU_logL_learnhypers']);
% cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_auxU_logL_learnhypers.eps']);
% system(cmd);

% figure; 
% plot(summaryMarg.LogL,'k'); 
% %xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% axis(ax);
% set(gca,'fontsize',fontsz);
% title('mGrad');
% set(gca, 'XTick', [3000 6000 9000]);
% print('-depsc2', '-r300', [outdir datasetName{1} '_marg_logL_learnhypers']);
% cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_marg_logL_learnhypers.eps']);
% system(cmd);

% figure; 
% plot(summaryEllipt.LogL,'k'); 
% xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% axis(ax);
% set(gca,'fontsize',fontsz);
% title('Ellipt-gibbs');
% %set(gca, 'XTick', [3000 6000 9000]);
% print('-depsc2', '-r300', [outdir datasetName '_ellipt_logL_learnhypers']);
% cmd = sprintf('epstopdf %s', [outdir datasetName '_ellipt_logL_learnhypers.eps']);
% system(cmd);

% figure; 
% plot(summarypCN.LogL,'k'); 
% xlabel('Sampling iteration','fontsize',fontsz);
% ylabel('Log likelihood','fontsize',fontsz);
% axis(ax);
% set(gca,'fontsize',fontsz);
% title('pCN');
% set(gca, 'XTick', [3000 6000 9000]);
% print('-depsc2', '-r300', [outdir datasetName{1} '_pCN_logL_learnhypers']);
% cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_pCN_logL_learnhypers.eps']);
% system(cmd);

% figure; 
% plot(summarypCNL.LogL,'k'); 
% xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% axis(ax);
% set(gca,'fontsize',fontsz);
% title('pCNL-gibbs');
% %set(gca, 'XTick', [3000 6000 9000]);
% print('-depsc2', '-r300', [outdir datasetName '_pCNL_logL_learnhypers']);
% cmd = sprintf('epstopdf %s', [outdir datasetName '_pCNL_logL_learnhypers.eps']);
% system(cmd);

%figure; 
%plot(summaryMALA.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
%axis(ax);
%set(gca,'fontsize',fontsz);
%title('pMALA');
%set(gca, 'XTick', [3000 6000 9000]);
%print('-depsc2', '-r300', [outdir datasetName{1} '_mala_logL_learnhypers']);
%cmd = sprintf('epstopdf %s', [outdir datasetName{1} '_mala_logL_learnhypers.eps']);
%system(cmd);


