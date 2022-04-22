function [results,stats] = getAllBurstAttr(results,LFP_seg,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',validStruct);
addRequired(p,'LFP_seg',@(x) validateattributes(x,{'cell','numeric'},{'vector'}));
addOptional(p,'butter_order',6,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addOptional(p,'ZS_threshold',2,@(x) validateattributes(x,{'numeric'},{'scalar'}));
% customized detection threshold
addParamValue(p,'threshold',[],validPositiveScalar);
% proportion of peak value at which burst duration stops
addParamValue(p,'stop_perc',0.25,@(x) validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
% window for FFT of bursts (>= 4 seconds)
addParamValue(p,'fft_length',4,validPositiveScalar);
addParamValue(p,'display',false,validTF);
addParamValue(p,'fit_ASamp',true,validTF);
addParamValue(p,'fit_ASamp_gamma',true,validTF);
addParamValue(p,'msgfcn',[]);
parse(p,results,LFP_seg,varargin{:});

fs = results.fs;
T_length = results.T_length;
valid_seg_id = results.valid_seg_id;

butter_order = p.Results.butter_order;
ZS_threshold = p.Results.ZS_threshold;
stop_perc = p.Results.stop_perc;
Nfft = 2^nextpow2(p.Results.fft_length*fs);

% Filter signal
% [bFilt,aFilt] = butter(butter_order,results.bump_f*2/fs);
[Z,P,K] = butter(butter_order,results.bump_f*2/fs);
[sos,g] = zp2sos(Z,P,K);    % using sos to avoid numeric error
nseg = numel(valid_seg_id);
LFP_ASamp = cell(nseg,1);
for i = 1:nseg
%     LFP_ASamp{i} = abs(hilbert(filtfilt(bFilt,aFilt,LFP_seg{valid_seg_id(i)})));
    LFP_ASamp{i} = abs(hilbert(filtfilt(sos,g,LFP_seg{valid_seg_id(i)})));
end
LFP_ASamp_all = cell2mat(LFP_ASamp);
mean_amp = mean(LFP_ASamp_all);
std_amp = std(LFP_ASamp_all);

% Filtered SNR
f = results.f;
fit_i = results.fit_i;
PSD_smoo = results.PSD_smoo;
PSD_bg = PSD_smoo;
PSD_fit = results.fit_a*results.f.^results.fit_b;
PSD_bg(fit_i) = PSD_fit(fit_i);
f_rsp = freqz(sos,f,fs)*g;	% using sos
m2_rsp = (f_rsp.*conj(f_rsp)).^2;	% square for filtfilt
pow_psd = @(x) trapz(f(2:end),m2_rsp(2:end).*x(2:end));
pow_bg_filt = pow_psd(PSD_bg);
pow_sig_filt = pow_psd(PSD_smoo)-pow_bg_filt;
SNR_filt = pow_sig_filt/pow_bg_filt;

% Calculate detection threshold
if isempty(p.Results.threshold)
    outliers = results.outliers;
    bump_i = results.bump_i;
    psd_cutoff = PSDcutoff(f,PSD_smoo,outliers, ...
        db2pow(results.nDB)*reshape(PSD_bg(outliers),size(outliers)));
%     psd_cutoff = PSDcutoff(f,PSD_smoo,outliers);
    pow_cutoff = trapz(f(bump_i(1):bump_i(2)),psd_cutoff(bump_i(1):bump_i(2)))-results.pow_bg;
    pow_cutoff_filt = pow_psd(psd_cutoff)-pow_bg_filt;
    prop = (pow_bg_filt+pow_cutoff_filt)/(pow_bg_filt+pow_sig_filt);
    LFP_ASamp_all = sort(LFP_ASamp_all);
    n = length(LFP_ASamp_all);
%     cumpow = cumsum(LFP_ASamp_all.^2);  % 2x cumulative energy
    cumpow = cumsum(LFP_ASamp_all.^2)./(1:n)';  % 2x cumulative power
    idx = max(find(cumpow>prop*cumpow(end),1)-1,1);
    threshold_auto = LFP_ASamp_all(idx);
    threshold = max(threshold_auto,mean_amp+ZS_threshold*std_amp);
else
    threshold = p.Results.threshold;
end

% Extract burst attributes
burst_seg = cell(nseg,1);
pks_seg_id = cell(nseg,1);
for i = 1:nseg
    burst_seg{i} = getBurstAttr(LFP_seg{valid_seg_id(i)},LFP_ASamp{i},fs,threshold,results.bump_f,stop_perc,Nfft);
    pks_seg_id{i} = int32(zeros(burst_seg{i}.npeak,1)+valid_seg_id(i));
end
burst_seg = cell2mat(burst_seg);
FieldsSum = {'npeak','ndetect','nmerged','nvalid'};
FieldsCat = {'peaks','tpks','AP','CN','BF'};
stats = burst_seg(1);
stats.pks_seg_id = cell2mat(pks_seg_id);
for i = 1:length(FieldsSum)
    stats.(FieldsSum{i}) = sum([burst_seg.(FieldsSum{i})]);
end
for i = 1:length(FieldsCat)
    stats.(FieldsCat{i}) = vertcat(burst_seg.(FieldsCat{i}));
end
ASAP_rate = stats.npeak/T_length;

if p.Results.fit_ASamp
    % AS amplitude distribution
    nbins = 60;
    n = length(LFP_ASamp_all);
    maxamp = min(max(LFP_ASamp_all),10*mean_amp);
    LFP_ASamp_all = LFP_ASamp_all(LFP_ASamp_all<=maxamp);
    as_amp = linspace(0,maxamp,nbins+1);
    as_amp = (as_amp(1:end-1)+as_amp(2:end))/2;
    hist_amp = hist(LFP_ASamp_all,as_amp);
    pdf_amp = hist_amp/n/binwidth(as_amp);
    
    % Fit gamma parameters for AS amplitude
    n = length(LFP_ASamp_all);
    if n>1e6
        rng(0);
        LFP_ASamp_all = LFP_ASamp_all(randi(n,[1e6,1]));
    end
    gam_amp_tot = gamfit(LFP_ASamp_all);
    if p.Results.fit_ASamp_gamma
        prec = 0.02;
        Grid = gamma_pdf_int_grid(gam_amp_tot(1),gam_amp_tot(2),prec,maxamp);
        mu_k = @(k) (pow_sig_filt*2/k/(k+1))^0.5;
        warning('off','all');
        gam_amp_sig = mle(LFP_ASamp_all,'pdf',@(x,k) ncxpdf_cond_gam(x,k,mu_k(k), ...
            pow_bg_filt^0.5,[],[],prec,Grid),'start',gam_amp_tot(1),'Lowerbound',0.1);
        warning('on','all');
        gam_amp_sig = [gam_amp_sig,mu_k(gam_amp_sig)];
    end
end

% Add results
varnames = {'butter_order','mean_amp','std_amp','stop_perc','ZS_threshold','threshold', ...
    'ASAP_rate','pow_sig_filt','pow_bg_filt','SNR_filt'};
if isempty(p.Results.threshold)
    varnames = [varnames,'pow_cutoff','pow_cutoff_filt','threshold_auto'];
end
if p.Results.fit_ASamp
    varnames = [varnames,'as_amp','hist_amp','pdf_amp','gam_amp_tot'];
    if p.Results.fit_ASamp_gamma
        varnames = [varnames,'gam_amp_sig'];
    end
end
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    if isempty(p.Results.threshold)
        dispfcn(num2str((threshold_auto-mean_amp)/std_amp,'Auto threshold = %.2f Z-score'));
    end
    dispfcn(num2str([threshold,(threshold-mean_amp)/std_amp],...
        ['Detection threshold = %.5g ',char(956),'V (%.2f Z-score)']));
    dispfcn(num2str(stats.npeak/T_length,'Total AS-AP rate = %.2f Hz'));
    dispfcn(num2str(stats.ndetect/T_length,'Detected AS-AP rate = %.2f Hz --- 100%%'));
    deducted = stats.ndetect-stats.nmerged;
    dispfcn(num2str([deducted/T_length,deducted/stats.ndetect*100],...
        'AS-AP rate deducted by merging = %.2f Hz --- %.2f%%'));
    deducted = stats.nmerged-stats.nvalid;
    dispfcn(num2str([deducted/T_length,deducted/stats.ndetect*100],...
        'AS-AP rate deducted by failed frequency estimation = %.2f Hz --- %.2f%%'));
    dispfcn(num2str([stats.nvalid/T_length,stats.nvalid/stats.ndetect*100],...
        'Valid detection rate = %.2f Hz --- %.2f%%'));
end

end