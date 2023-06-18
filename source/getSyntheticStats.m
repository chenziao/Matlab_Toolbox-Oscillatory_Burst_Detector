function [SynStats,SynData] = getSyntheticStats(SynParam,SynData)
%% Synthetic Data
SynData.CompositeTrace = SynData.BackgroundTrace+SynData.SignalTrace;
syn_len = SynParam.syn_len;
NT = SynParam.NT;
unit = SynParam.rdat.unit;

%% Primary Process
LFP_filt = zeros(NT,3);
LFP_amp = zeros(NT,3);
peaks = cell(1,3);  tpks = cell(1,3);
npk = zeros(1,3);   AP_rates = zeros(1,3);
fieldname = {'BackgroundTrace','SignalTrace','CompositeTrace'};
for i = 1:3
    LFP_filt(:,i) = filtfilt(SynParam.sos,SynParam.g,SynData.(fieldname{i}));   % using sos
    LFP_amp(:,i) = abs(hilbert(LFP_filt(:,i)));
    [peaks{i},tpks{i}] = Findpeaks(LFP_amp(:,i));
    npk(i) = length(peaks{i});
    AP_rates(i) = npk(i)/syn_len;
end

%% Statistics
sigma_syn = std(LFP_filt,1,1);
max_amp = max(LFP_amp,[],1);
amp_mean = mean(LFP_amp,1);
amp_sigma = std(LFP_amp,1,1);
amp_skew = skewness(LFP_amp(:,2),1,1);
amp_kurt = kurtosis(LFP_amp(:,2),1,1);
gam_amp_syn = cell2mat(arrayfun(@(i) gamfit(LFP_amp(:,i)),1:3,'UniformOutput',false)');
gam_AP_syn = cell2mat(cellfun(@gamfit,peaks,'UniformOutput',false)');
gam_amp_tot = SynParam.rdat.gam_amp_tot;
gam_AP_tot = SynParam.rdat.gam_AP;

%% Output
varnames = {'LFP_amp','unit','syn_len','NT','peaks','tpks','npk','AP_rates', ...
    'sigma_syn','max_amp','amp_mean','amp_sigma','amp_skew','amp_kurt', ...
    'gam_amp_syn','gam_AP_syn','gam_amp_tot','gam_AP_tot'};
SynStats = struct();
for i = 1:numel(varnames),	SynStats.(varnames{i})=eval(varnames{i});	end

end
