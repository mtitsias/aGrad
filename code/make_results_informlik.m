outdir = '../diagrams/';
fontsz = 26;
addpath ../results/; 

load ../data/regressinformlik.mat; 

%load regressionInformLikeToy_AuxU.mat; 
%load regressionInformLikeToy_AuxZ.mat;
%load regressionInformLikeToy_Marg.mat;
%load regressionInformLikeToy_MALA.mat;
%load regressionInformLikeToy_Ellipt.mat;
%load regressionInformLikeToy_pCN.mat;
%load regressionInformLikeToy_pCNL.mat;

datasetName = 'RegressInformLikel';

Repeats = 10; 

for j=1:3
     
figure;
plot(XX{j},YY{j},'ko','markersize',3);
xlabel('Input','fontsize',fontsz);
if j == 1
ylabel('Response','fontsize',fontsz);
end
set(gca,'fontsize',fontsz);
name = [outdir datasetName '_data' num2str(j)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);   

for rep=1:Repeats

load(['regressionInformLikeToy_AuxU_repeat' num2str(rep) '.mat']); 
load(['regressionInformLikeToy_AuxZ_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_Marg_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_MALA_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_Ellipt_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_pCN_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_pCNL_repeat' num2str(rep) '.mat']);  
    
ESSminAuxZ(rep) = min(summaryAuxZ{j}.eff_F); 
ESSmedianAuxZ(rep) = median(summaryAuxZ{j}.eff_F); 
ESSmaxAuxZ(rep) = max(summaryAuxZ{j}.eff_F); 
ESSlogLAuxZ(rep) = summaryAuxZ{j}.eff_LogL; 
TrainTimeAuxZ(rep) = summaryAuxZ{j}.elapsed;
deltaAuxZ(rep) = summaryAuxZ{j}.delta; 

ESSminAuxU(rep) = min(summaryAuxU{j}.eff_F); 
ESSmedianAuxU(rep) = median(summaryAuxU{j}.eff_F); 
ESSmaxAuxU(rep) = max(summaryAuxU{j}.eff_F); 
ESSlogLAuxU(rep) = summaryAuxU{j}.eff_LogL; 
TrainTimeAuxU(rep) = summaryAuxU{j}.elapsed;
deltaAuxU(rep) = summaryAuxU{j}.delta; 

ESSminMarg(rep) = min(summaryMarg{j}.eff_F); 
ESSmedianMarg(rep) = median(summaryMarg{j}.eff_F); 
ESSmaxMarg(rep) = max(summaryMarg{j}.eff_F); 
ESSlogLMarg(rep) = summaryMarg{j}.eff_LogL; 
TrainTimeMarg(rep) = summaryMarg{j}.elapsed;
deltaMarg(rep) = summaryMarg{j}.delta;

ESSminMALA(rep) = min(summaryMALA{j}.eff_F); 
ESSmedianMALA(rep) = median(summaryMALA{j}.eff_F); 
ESSmaxMALA(rep) = max(summaryMALA{j}.eff_F); 
ESSlogLMALA(rep) = summaryMALA{j}.eff_LogL; 
TrainTimeMALA(rep) = summaryMALA{j}.elapsed;
deltaMALA(rep) = summaryMALA{j}.delta;

ESSminEllipt(rep) = min(summaryEllipt{j}.eff_F); 
ESSmedianEllipt(rep) = median(summaryEllipt{j}.eff_F); 
ESSmaxEllipt(rep) = max(summaryEllipt{j}.eff_F); 
ESSlogLEllipt(rep) = summaryEllipt{j}.eff_LogL; 
TrainTimeEllipt(rep) = summaryEllipt{j}.elapsed;

ESSminpCN(rep) = min(summarypCN{j}.eff_F); 
ESSmedianpCN(rep) = median(summarypCN{j}.eff_F); 
ESSmaxpCN(rep) = max(summarypCN{j}.eff_F); 
ESSlogLpCN(rep) = summarypCN{j}.eff_LogL; 
TrainTimepCN(rep) = summarypCN{j}.elapsed;
deltapCN(rep) = summarypCN{j}.delta;

ESSminpCNL(rep) = min(summarypCNL{j}.eff_F); 
ESSmedianpCNL(rep) = median(summarypCNL{j}.eff_F); 
ESSmaxpCNL(rep) = max(summarypCNL{j}.eff_F); 
ESSlogLpCNL(rep) = summarypCNL{j}.eff_LogL; 
TrainTimepCNL(rep) = summarypCNL{j}.elapsed;
deltapCNL(rep) = summarypCNL{j}.delta;

if rep == 1

figure; 
plot(summaryAuxZ{j}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
if j==1
ylabel('Log-likelihood','fontsize',fontsz);
end
set(gca,'fontsize',fontsz);
if (j==1)
axis([0 15000 -1550 -1400]);
elseif (j==2)
axis([0 15000 -360 -240]); 
elseif (j==3) 
axis([0 15000 800 900]);
elseif (j==4)
axis([0 15000 1900 2060]); 
end
box on;
title('aGrad-z');
name = [outdir datasetName '_auxZ_logL' num2str(j)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

figure; 
plot(summaryAuxU{j}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
if j==1
ylabel('Log-likelihood','fontsize',fontsz);
end
set(gca,'fontsize',fontsz);
if (j==1)
axis([0 15000 -1550 -1400]);
elseif (j==2)
axis([0 15000 -360 -240]); 
elseif (j==3) 
axis([0 15000 800 900]);
elseif (j==4)
axis([0 15000 1900 2060]); 
end
box on; 
title('aGrad-u');
name = [outdir datasetName '_auxU_logL' num2str(j)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

figure; 
hold on; 
plot(summaryMarg{j}.LogL,'k'); 
%xlabel('Sampling iteration','fontsize',fontsz);
if j==1
ylabel('Log-likelihood','fontsize',fontsz);
end
set(gca,'fontsize',fontsz);
if (j==1)
axis([0 15000 -1550 -1400]);
elseif (j==2)
axis([0 15000 -360 -240]); 
elseif (j==3) 
axis([0 15000 800 900]);
elseif (j==4)
axis([0 15000 1900 2060]); 
end
box on;
title('mGrad');
name = [outdir datasetName '_marg_logL' num2str(j)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

figure; 
hold on;
plot(summaryEllipt{j}.LogL,'b'); 
xlabel('Sampling iteration','fontsize',fontsz);
if j==1
ylabel('Log-likelihood','fontsize',fontsz);
end
set(gca,'fontsize',fontsz);
if (j==1)
axis([0 15000 -1550 -1400]);
elseif (j==2)
axis([0 15000 -360 -240]); 
elseif (j==3) 
axis([0 15000 800 900]);
elseif (j==4)
axis([0 15000 1900 2060]); 
end
%title('Ellipt');
%name = [outdir datasetName '_ellipt_logL' num2str(j)];
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);

%figure; 
plot(summarypCN{j}.LogL,'r'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
%set(gca,'fontsize',fontsz);
%if (j==1)
%axis([0 10000 -1550 -1400]);
%elseif (j==2)
%axis([0 10000 -340 -240]); 
%elseif (j==3) 
%axis([0 10000 0 900]);
%elseif (j==4)
%axis([0 10000 0 2060]); 
%end
%title('pCN');
%name = [outdir datasetName '_pCN_logL' num2str(j)];
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);

%figure; 
plot(summarypCNL{j}.LogL,'g'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
%set(gca,'fontsize',fontsz);
%if (j==1)
%axis([0 10000 -1550 -1400]);
%elseif (j==2)
%axis([0 10000 -340 -240]); 
%elseif (j==3) 
%axis([0 10000 0 900]);
%elseif (j==4)
%axis([0 10000 0 2060]); 
%end
%title('pCNL');
%name = [outdir datasetName '_pCNL_logL' num2str(j)];
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);


%figure; 
plot(summaryMALA{j}.LogL,'m'); 
%xlabel('Sampling iteration','fontsize',fontsz);
%ylabel('Log likelihood','fontsize',fontsz);
%set(gca,'fontsize',fontsz);
%if (j==1)
%axis([0 10000 -1550 -1400]);
%elseif (j==2)
%axis([0 10000 -340 -240]); 
%elseif (j==3) 
%axis([0 10000 0 900]);
%elseif (j==4)
%axis([0 10000 0 2060]); 
%end
%title('MALA');
%name = [outdir datasetName '_mala_logL' num2str(j)];
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);
box on; 
%title('MALA');
if j==1
legend('Ellipt','pCN','pCNL', 'pMALA', 'Location','SouthEast');
end

name = [outdir datasetName '_rest_logL' num2str(j)];
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

end


end

fid = fopen([outdir datasetName '_table' num2str(j) '.txt'],'w'); 
fprintf(fid,'Method &  Time(s) & Step size $\\delta$  &  ESS (Min, Med, Max)  & Min ESS/s (1 st.d.)\\\\ \n'); 
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