outdir = '../diagrams/';
fontsz = 26;
addpath ../results/; 

 
datasetName = 'GaussianCox';

Repeats = 10;

for i=1:2

for rep=1:Repeats

load(['logGaussianCoxGirolami_AuxZ_repeat' num2str(rep) '.mat']);

ESSminAuxZ(rep) = min(summaryAuxZ{i}.eff_F); 
ESSmedianAuxZ(rep) = median(summaryAuxZ{i}.eff_F); 
ESSmaxAuxZ(rep) = max(summaryAuxZ{i}.eff_F); 
ESSlogLAuxZ(rep) = summaryAuxZ{i}.eff_LogL; 
TrainTimeAuxZ(rep) = summaryAuxZ{i}.elapsed;
deltaAuxZ(rep) = summaryAuxZ{i}.delta; 

load(['logGaussianCoxGirolami_AuxU_repeat' num2str(rep) '.mat']);
ESSminAuxU(rep) = min(summaryAuxU{i}.eff_F); 
ESSmedianAuxU(rep) = median(summaryAuxU{i}.eff_F); 
ESSmaxAuxU(rep) = max(summaryAuxU{i}.eff_F); 
ESSlogLAuxU(rep) = summaryAuxU{i}.eff_LogL; 
TrainTimeAuxU(rep) = summaryAuxU{i}.elapsed;
deltaAuxU(rep) = summaryAuxU{i}.delta; 

load(['logGaussianCoxGirolami_Marg_repeat' num2str(rep) '.mat']);
ESSminMarg(rep) = min(summaryMarg{i}.eff_F); 
ESSmedianMarg(rep) = median(summaryMarg{i}.eff_F); 
ESSmaxMarg(rep) = max(summaryMarg{i}.eff_F); 
ESSlogLMarg(rep) = summaryMarg{i}.eff_LogL; 
TrainTimeMarg(rep) = summaryMarg{i}.elapsed;
deltaMarg(rep) = summaryMarg{i}.delta;

load(['logGaussianCoxGirolami_mala_repeat' num2str(rep) '.mat']);
ESSminMALA(rep) = min(summaryMALA{i}.eff_F); 
ESSmedianMALA(rep) = median(summaryMALA{i}.eff_F); 
ESSmaxMALA(rep) = max(summaryMALA{i}.eff_F); 
ESSlogLMALA(rep) = summaryMALA{i}.eff_LogL; 
TrainTimeMALA(rep) = summaryMALA{i}.elapsed;
deltaMALA(rep) = summaryMALA{i}.delta;

load(['logGaussianCoxGirolami_Ellipt_repeat' num2str(rep) '.mat']);
ESSminEllipt(rep) = min(summaryEllipt{i}.eff_F); 
ESSmedianEllipt(rep) = median(summaryEllipt{i}.eff_F); 
ESSmaxEllipt(rep) = max(summaryEllipt{i}.eff_F); 
ESSlogLEllipt(rep) = summaryEllipt{i}.eff_LogL; 
TrainTimeEllipt(rep) = summaryEllipt{i}.elapsed;

load(['logGaussianCoxGirolami_pCN_repeat' num2str(rep) '.mat']);
ESSminpCN(rep) = min(summarypCN{i}.eff_F); 
ESSmedianpCN(rep) = median(summarypCN{i}.eff_F); 
ESSmaxpCN(rep) = max(summarypCN{i}.eff_F); 
ESSlogLpCN(rep) = summarypCN{i}.eff_LogL; 
TrainTimepCN(rep) = summarypCN{i}.elapsed;
deltapCN(rep) = summarypCN{i}.delta;

load(['logGaussianCoxGirolami_pCNL_repeat' num2str(rep) '.mat']);
ESSminpCNL(rep) = min(summarypCNL{i}.eff_F); 
ESSmedianpCNL(rep) = median(summarypCNL{i}.eff_F); 
ESSmaxpCNL(rep) = max(summarypCNL{i}.eff_F); 
ESSlogLpCNL(rep) = summarypCNL{i}.eff_LogL; 
TrainTimepCNL(rep) = summarypCNL{i}.elapsed;
deltapCNL(rep) = summarypCNL{i}.delta;

% load logGaussianCoxGirolami_mMALA_RHMC.mat;
% ESSminRHMC = min(summaryRHMC{i}.eff_F); 
% ESSmedianRHMC = median(summaryRHMC{i}.eff_F); 
% ESSmaxRHMC = max(summaryRHMC{i}.eff_F); 
% ESSlogLRHMC = summaryRHMC{i}.eff_LogL; 
% TrainTimeRHMC = summaryRHMC{i}.elapsed;
% 
% ESSminmMALA = min(summarymMALA{i}.eff_F); 
% ESSmedianmMALA = median(summarymMALA{i}.eff_F); 
% ESSmaxmMALA = max(summarymMALA{i}.eff_F); 
% ESSlogLmMALA = summarymMALA{i}.eff_LogL; 
% TrainTimemMALA = summarymMALA{i}.elapsed;

if rep == 1

figure; 
plot(summaryAuxZ{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('aGrad-z');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_auxZ_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


figure; 
plot(summaryAuxU{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('aGrad-u');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_auxU_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


figure; 
plot(summaryMarg{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
axis([1 7000 2150 2250]);
set(gca,'fontsize',fontsz);
title('mGrad');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_marg_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


figure; 
plot(summaryEllipt{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('Ellipt');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_ellipt_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


figure; 
plot(summarypCN{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('pCN');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_pCN_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


figure; 
plot(summarypCNL{i}.LogL,'k'); 
xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('pCNL');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_pCNL_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

figure; 
plot(summaryMALA{i}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
ylabel('Log-likelihood','fontsize',fontsz);
set(gca,'fontsize',fontsz);
axis([1 7000 2150 2250]);
title('pMALA');
set(gca, 'XTick', [2000 4000 6000]);
name = [outdir datasetName '_mala_logL' num2str(i)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);


% figure; 
% plot(summarymMALA{i}.LogL,'k'); 
% xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% set(gca,'fontsize',fontsz);
% axis([1 7000 2150 2250]);
% title('mMALA');
% set(gca, 'XTick', [2000 4000 6000]);
% name = [outdir datasetName '_mmala_logL'];
% print('-depsc2', '-r300', name);
% cmd = sprintf('epstopdf %s', [name '.eps']);
% system(cmd);
% figure; 
% plot(summaryRHMC{i}.LogL,'k'); 
% xlabel('Sampling iteration','fontsize',fontsz);
% %ylabel('Log likelihood','fontsize',fontsz);
% set(gca,'fontsize',fontsz);
% axis([1 7000 2150 2250]);
% title('RHMC');
% set(gca, 'XTick', [2000 4000 6000]);
% name = [outdir datasetName '_RHMC_logL'];
% print('-depsc2', '-r300', name);
% cmd = sprintf('epstopdf %s', [name '.eps']);
% system(cmd);

end

end

fid = fopen([outdir datasetName '_table' num2str(i) '.txt'],'w'); 
fprintf(fid,'Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.) \\\\ \n'); 
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


end



figure;  
imagesc(reshape(summaryAuxZ{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxZ_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxZ_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summaryAuxZ{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxZ_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxZ_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summaryAuxZ{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxZ_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxZ_expfield.eps']);
system(cmd);


figure;  
imagesc(reshape(summaryAuxU{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxU_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxU_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summaryAuxU{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxU_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxU_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summaryAuxU{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_auxU_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_auxU_expfield.eps']);
system(cmd);

figure;  
imagesc(reshape(summaryMarg{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_marg_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_marg_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summaryMarg{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_marg_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_marg_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summaryMarg{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_marg_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_marg_expfield.eps']);
system(cmd);

figure;  
imagesc(reshape(summaryEllipt{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_ellipt_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_ellipt_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summaryEllipt{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_ellipt_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_ellipt_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summaryEllipt{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_ellipt_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_ellipt_expfield.eps']);
system(cmd);


figure;  
imagesc(reshape(summarypCN{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCN_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCN_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summarypCN{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCN_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCN_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summarypCN{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCN_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCN_expfield.eps']);
system(cmd);


figure;  
imagesc(reshape(summarypCNL{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCNL_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCNL_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summarypCNL{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCNL_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCNL_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summarypCNL{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_pCNL_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_pCNL_expfield.eps']);
system(cmd);

figure;  
imagesc(reshape(summaryMALA{1}.meanF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_mala_meanfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mala_meanfield.eps']);
system(cmd);
figure;  
imagesc(reshape(summaryMALA{1}.varF,64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_mala_varfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mala_varfield.eps']);
system(cmd);
figure;  
imagesc(reshape(exp(summaryMALA{1}.meanF),64,64)); 
axis off;
print('-depsc2', '-r300', [outdir 'GaussianCox_mala_expfield']);
cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mala_expfield.eps']);
system(cmd);


% figure;  
% imagesc(reshape(summarymMALA{1}.meanF,64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_mmala_meanfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mmala_meanfield.eps']);
% system(cmd);
% figure;  
% imagesc(reshape(summarymMALA{1}.varF,64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_mmala_varfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mmala_varfield.eps']);
% system(cmd);
% figure;  
% imagesc(reshape(exp(summarymMALA{1}.meanF),64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_mmala_expfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_mmala_expfield.eps']);
% system(cmd);
% figure;  
% imagesc(reshape(summaryRHMC{1}.meanF,64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_RHMC_meanfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_RHMC_meanfield.eps']);
% system(cmd);
% figure;  
% imagesc(reshape(summaryRHMC{1}.varF,64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_RHMC_varfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_RHMC_varfield.eps']);
% system(cmd);
% figure;  
% imagesc(reshape(exp(summaryRHMC{1}.meanF),64,64)); 
% axis off;
% print('-depsc2', '-r300', [outdir 'GaussianCox_RHMC_expfield']);
% cmd = sprintf('epstopdf %s', [outdir 'GaussianCox_RHMC_expfield.eps']);
% system(cmd);



