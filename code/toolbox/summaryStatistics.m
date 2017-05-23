function summary = summaryStatistics(samples)
%
%

% compute effective sample sizes
summary.eff_F = zeros(1, size(samples.F,2));

%for j=1:floor(size(samples.F,2)/8)
%    set = (8*(j-1)+1):(j*8);
%    summary.eff_F(set) = effective_size_rcoda(samples.F(:,set));
%end
%if (8*j+1) < size(samples.F,2)
%   set = (8*j+1):size(samples.F,2);
%   summary.eff_F(set) = effective_size_rcoda(samples.F(:,set));
%end

for j=1:size(samples.F,2)
    summary.eff_F(j) = mcmc_ess(samples.F(:,j));
end
%if (8*j+1) < size(samples.F,2)
%   set = (8*j+1):size(samples.F,2);
%   summary.eff_F(set) = effective_size_rcoda(samples.F(:,set));
%end

%summary.eff_LogL = mcmc_ess(samples.LogL(:));

% compute mean and variance 
summary.meanF = mean(samples.F);
summary.varF = var(samples.F);

% keep the log likelihoods 
summary.LogL = samples.LogL;

if isfield(samples,'num_calls')
   summary.num_calls = samples.num_calls;
end

