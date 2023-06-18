function BoundStats = getPreselectedBounds(SynParam,SynStats,tL,tU,auto_thetaU)
% tL,tU - Lower/Upper bounds in multiple of sigma of background trace.
% auto_thetaU - to keep track of value in GUI, not for customized usage.
if nargin<3 || isempty(tL)
    tL = 2^0.5;
end
if nargin<4
    tU = [];
end
if nargin<5
    auto_thetaU = [];
end
LFP_amp = SynStats.LFP_amp(:,2);
tpks = SynStats.tpks;
syn_len = SynStats.syn_len;
NT = SynStats.NT;
AP_rates = SynStats.AP_rates;
sigma = SynStats.sigma_syn(1);

%% Calculate Lower/Upper bounds
if isempty(auto_thetaU)
    % Power in filtered band
    [psd_bg,f,~,~,psd_tot] = SynParam.psd_background;
    f_rsp = freqz(SynParam.sos,f,SynParam.fs)*SynParam.g;	% using sos
    m2_rsp = (f_rsp.*conj(f_rsp)).^2;	% square for filtfilt
    pow_psd = @(x) trapz(f(2:end),m2_rsp(2:end).*x(2:end));
    
    % Calculate power below cutoff dB threshold
    outliers = SynParam.rdat.outliers;
    outliers = round(interp1(f,1:numel(f), ...
        reshape(SynParam.rdat.f(outliers),size(outliers))));
    psd_cutoff = PSDcutoff(psd_tot,outliers,db2pow(SynParam.rdat.tDB)*psd_bg);
    pow_T = pow_psd(psd_cutoff-psd_bg);
    prop = pow_T/pow_psd(psd_tot-psd_bg);
    
    % Find cutoff amplitude
    LFP_amp_sort = sort(LFP_amp);
    % cumpow = cumsum(LFP_amp_sort.^2)/SynStats.NT/2; % cumulative energy
    cumpow = cumsum(LFP_amp_sort.^2)./(1:NT)'/2;	% cumulative power
    idx = max(find(cumpow>prop*cumpow(end),1)-1,1); % use relative power
    % idx = max(find(cumpow>pow_T,1)-1,1);    % use absolute power
    auto_thetaU = (LFP_amp_sort(idx)+LFP_amp_sort(idx+1))/2;
end

pre_bound = [0,0];
% Lower bound
pre_bound(1) = sigma*tL;
% Upper bound
if isempty(tU)
    pre_bound(2) = auto_thetaU;
else
    pre_bound(2) = sigma*tU;
end
pre_bound(2) = max(pre_bound);	% Restrict Lower<=Upper
nsigma_bound = pre_bound/sigma;
ZS_bound = (pre_bound-SynStats.amp_mean(1))/SynStats.amp_sigma(1);

%% Burst rate
% Proportion of each duration
F_idx = LFP_amp<pre_bound(1);
T_idx = LFP_amp>pre_bound(2);
dur_prop = [sum(F_idx),0,sum(T_idx)]/NT;
dur_prop(2) = 1-sum(dur_prop);
% AP rate in composite trace
APrate_fit = [sum(F_idx(tpks{3})),0,sum(T_idx(tpks{3}))]/syn_len;
APrate_fit(2) = AP_rates(3)-sum(APrate_fit);
% AP rate in signal trace
APrate_sig_fit = [sum(F_idx(tpks{2})),0,sum(T_idx(tpks{2}))]/syn_len;
APrate_sig_fit(2) = AP_rates(2)-sum(APrate_sig_fit);

%% Output
varnames = {'pre_bound','nsigma_bound','ZS_bound','auto_thetaU', ...
    'F_idx','T_idx','dur_prop','APrate_fit','APrate_sig_fit'};
BoundStats = struct();
for i = 1:numel(varnames),	BoundStats.(varnames{i})=eval(varnames{i});	end

end
