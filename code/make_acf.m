outdir = '../diagrams/';
fontsz = 26;
addpath ../results/; 

load ../data/regressinformlik.mat; 
datasetName = 'RegressInformLikel';

numLags = 4000; 
datasetName = 'RegressInformLikel';
Repeats = 10; 

for j=1:3
    
acfAuxU = zeros(numLags+1, Repeats);
acfAuxZ = zeros(numLags+1, Repeats);
acfMarg = zeros(numLags+1, Repeats);
acfEllipt = zeros(numLags+1, Repeats);
acfMALA = zeros(numLags+1, Repeats);
acfpCN = zeros(numLags+1, Repeats);
acfpCNL = zeros(numLags+1, Repeats);

for rep=1:Repeats

load(['regressionInformLikeToy_AuxU_repeat' num2str(rep) '.mat']); 
load(['regressionInformLikeToy_AuxZ_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_Marg_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_MALA_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_Ellipt_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_pCN_repeat' num2str(rep) '.mat']);
load(['regressionInformLikeToy_pCNL_repeat' num2str(rep) '.mat']);  

[acfAuxU(:,rep), lags] = autocorr(summaryAuxU{j}.LogL(10001:end),numLags);
cpuTimeAuxU(rep) = summaryAuxU{j}.elapsed;
[acfAuxZ(:,rep), lags] = autocorr(summaryAuxZ{j}.LogL(10001:end),numLags);
cpuTimeAuxZ(rep) = summaryAuxZ{j}.elapsed;
[acfMarg(:,rep), lags] = autocorr(summaryMarg{j}.LogL(10001:end),numLags);
cpuTimeMarg(rep) = summaryMarg{j}.elapsed;


[acfEllipt(:,rep), lags] = autocorr(summaryEllipt{j}.LogL(end-5000:end),numLags);
cpuTimeEllipt(rep) = summaryEllipt{j}.elapsed;
[acfMALA(:,rep), lags] = autocorr(summaryMALA{j}.LogL(end-5000:end),numLags);
cpuTimeMALA(rep) = summaryMALA{j}.elapsed;
[acfpCN(:,rep), lags] = autocorr(summarypCN{j}.LogL(end-5000:end),numLags);
cpuTimepCN(rep) = summarypCN{j}.elapsed;

[acfpCNL(:,rep),lags] = autocorr(summarypCNL{j}.LogL(end-5000:end),numLags);
cpuTimepCNL(rep) = summarypCNL{j}.elapsed;

end

totalIters = 15000;
numLags = 100;
plotnumLags = 20;

figure; 
hold on;
plot(lags(1:numLags)*mean(cpuTimeAuxU)/totalIters, mean(acfAuxU(1:numLags,:),2),'ko-');
plot(lags(1:numLags)*mean(cpuTimeAuxZ)/totalIters, mean(acfAuxZ(1:numLags,:),2),'k*-');
plot(lags(1:numLags)*mean(cpuTimeMarg)/totalIters, mean(acfMarg(1:numLags,:),2),'kd-');
plot(lags(1:numLags)*mean(cpuTimeEllipt)/totalIters, mean(acfEllipt(1:numLags,:),2),'bo-');
plot(lags(1:numLags)*mean(cpuTimepCN)/totalIters, mean(acfpCN(1:numLags,:),2),'ro-');
plot(lags(1:numLags)*mean(cpuTimepCNL)/totalIters, mean(acfpCNL(1:numLags,:),2),'go-');
plot(lags(1:numLags)*mean(cpuTimeMALA)/totalIters, mean(acfMALA(1:numLags,:),2),'mo-');
axis([0 plotnumLags*(mean(cpuTimeAuxZ)/totalIters) 0 1]);
grid('on');
box on; 
if j==1
legend('aGrad-u','aGrad-z','mGrad','Ellipt','pCN','pCNL','pMALA', 'Location', 'best');%...
end

name = [outdir datasetName '_acf' num2str(j)];
xlabel('CPU time')
ylabel('Autocorrelation');
set(gca,'fontsize',fontsz);
print('-depsc2', '-r300', name);
cmd = sprintf('epstopdf %s', [name '.eps']);
system(cmd);

%figure; 
%hold on;
%plot(lags(1:numLags)*mean(cpuTimeEllipt)/totalIters, mean(acfEllipt(1:numLags,:),2),'bo-');
%plot(lags(1:numLags)*mean(cpuTimepCN)/totalIters, mean(acfpCN(1:numLags,:),2),'ro-');
%plot(lags(1:numLags)*mean(cpuTimepCNL)/totalIters, mean(acfpCNL(1:numLags,:),2),'go-');
%plot(lags(1:numLags)*mean(cpuTimeMALA)/totalIters, mean(acfMALA(1:numLags,:),2),'mo-');
%axis([0 plotnumLags*(mean(cpuTimeAuxZ)/totalIters) 0 1]);
%grid('on');
%box on; 
%legend('Ellipt','pCN','pCNL', 'pMALA', 'Location','SouthWest');
%name = [outdir datasetName '_rest_acf' num2str(j)];
%xlabel('CPU time')
%ylabel('Autocorrelation');
%set(gca,'fontsize',fontsz);
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);


%titleStr =  'aGrad-u';
%name = [outdir datasetName '_auxU_acf' num2str(j)];
%CreatACFplot(lags(1:100), mean(cpuTimeAuxU), mean(acfAuxU(1:100,:),2), numMA, numLags, bounds, name, titleStr, fontsz);
%
%titleStr =  'aGrad-z';
%name = [outdir datasetName '_auxZ_acf' num2str(j)];
%CreatACFplot(lags(1:100), mean(cpuTimeAuxZ), mean(acfAuxZ(1:100,:),2), numMA, numLags, bounds, name, titleStr, fontsz);
%
%titleStr =  'mGrad';
%name = [outdir datasetName '_marg_acf' num2str(j)];
%CreatACFplot(lags(1:100), mean(cpuTimeMarg), mean(acfMarg(1:100,:),2), numMA, numLags, bounds, name, titleStr, fontsz);

end


%function e = CreatACFplot(lags, cpuTime, acf, numMA, numLags, bounds, name, titleStr, fontsz)
%
%lineHandles = stem(lags*cpuTime, acf,'filled','k-o');
%axis([0 99 0 1]);
%set(lineHandles(1),'MarkerSize',4)
%grid('on')
%xlabel('Iteration scaled by total CPU time')
%ylabel('Autocorrelation')
%title('Sample Autocorrelation Function')
%hold('on');
%plot([numMA+0.5 numMA+0.5; numLags numLags],[bounds([1 1]) bounds([2 2])],'-b');
%plot([0 numLags],[0 0],'-k');
%box on;
%title(titleStr);
%set(gca,'fontsize',fontsz);
%print('-depsc2', '-r300', name);
%cmd = sprintf('epstopdf %s', [name '.eps']);
%system(cmd);
%
%hold('off')
%a = axis;
%axis([a(1:3) 1]);