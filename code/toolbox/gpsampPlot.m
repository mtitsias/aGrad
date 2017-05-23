function gpsampPlot(model, samples)
%
% 



% plot likelihood hyperparameters 
if isfield(samples,'likLogtheta')
   %
   for i=1:model.Likelihood.nParams; 
       figure;
       hist(exp(samples.likLogtheta(:,i)),50); 
       titlestring = 'likelihood hyperp: ';
       titlestring = [titlestring, num2str(i)]; 
       title(titlestring,'fontsize', 20);
   end
   %
end


% plot kernel hyperparameters 
if isfield(samples, 'kernLogtheta')
   %
   for i=1:model.GP.nParams; 
       figure;
       hist(exp(samples.kernLogtheta(:,i)),50); 
       titlestring = 'kernel hyperp: ';
       titlestring = [titlestring, num2str(i)]; 
       title(titlestring,'fontsize', 20);
   end
   %
end


% plot GP function 
if model.D == 1
   %
   %[sortX inX] = sort(model.X);
  
   if strcmp(model.Likelihood.type,'Poisson')
      mu = mean(exp(samples.F))';
      stds = sqrt(var(exp(samples.F)))';
   else
      mu = mean(samples.F)';
      stds = sqrt(var(samples.F))';
   end
   
   %stdsWithnoise = sqrt(var(samples.F))';
    
   figure
   hold on;
   fillColor = [0.7 0.7 0.7];
   %fillColor = [0.8 0.8 0.8];  % for the paper
   fill([model.X; model.X(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
   plot(model.X, mu,'b','lineWidth',3);
     
   %if strcmp(model.Likelihood.type,'Gaussian')
   plot(model.X, model.y, '+k', 'lineWidth', 1);
   %end
   %
end

