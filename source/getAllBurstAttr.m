function [results,stats] = getAllBurstAttr(results,LFP_seg,varargin)
p = inputParser;
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',@(x) validateattributes(x,{'struct'},{'scalar'}));
addRequired(p,'LFP_seg',@(x) validateattributes(x,{'cell','numeric'},{'vector'}));
% butterworth filter order
addOptional(p,'butter_order',6,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
% amplitude Z-score threshold for burst detection
addOptional(p,'ZS_threshold',2,@(x) validateattributes(x,{'numeric'},{'scalar'}));
% automatic threshold (not recommended)
addParameter(p,'auto_threshold',false,validTF);
% fixed amplitude threshold
addParameter(p,'threshold',[],@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
% proportion of peak value at which burst duration stops
addParameter(p,'stop_perc',0.25,@(x) validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
% resolution for burst frequency (<= 0.25 Hz)
addParameter(p,'freq_resolution',0.25,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'display',false,validTF);
% fit amplitude distribution
addParameter(p,'fit_amp',true,validTF);
% 0: Do not; 1: fit signal component amplitude distribution; 2: considering background power scale.
addParameter(p,'fit_amp_gamma',true,validTF);
% unit of LFP
addParameter(p,'unit','μV',@(x) validateattributes(x,{'char','string'},{'scalartext','nonempty'}));
addParameter(p,'msgfcn',[]);
parse(p,results,LFP_seg,varargin{:});

fs = results.fs;
T_length = results.T_length;
valid_seg_id = results.valid_seg_id;

butter_order = p.Results.butter_order;
stop_perc = p.Results.stop_perc;
freq_resolution = p.Results.freq_resolution;
Nfft = max(2^nextpow2(fs/freq_resolution),32);  % at least 32
unit = char(p.Results.unit);

%% Filter signal
% [bFilt,aFilt] = butter(butter_order,results.sig_range_f*2/fs);
[Z,P,K] = butter(butter_order,results.sig_range_f*2/fs);
[sos,g] = zp2sos(Z,P,K);    % using sos for numerical stability
nseg = numel(valid_seg_id);
LFP_AS = cell(nseg,1);  % analytic signal
LFP_amp = LFP_AS;
for i = 1:nseg
%     LFP_AS{i} = hilbert(filtfilt(bFilt,aFilt,LFP_seg{valid_seg_id(i)}));
    LFP_AS{i} = hilbert(filtfilt(sos,g,LFP_seg{valid_seg_id(i)}));
    LFP_amp{i} = abs(LFP_AS{i});
end
LFP_amp_all = cell2mat(LFP_amp);
mean_amp = mean(LFP_amp_all);
std_amp = std(LFP_amp_all);

%% Filtered SNR
f = results.f;
bg_range_i = results.bg_range_i;
PSD_smoo = results.PSD_smoo;
PSD_bg = PSD_smoo;
PSD_fit = getBackgroundFit(results);
PSD_bg(bg_range_i) = PSD_fit(bg_range_i);
f_rsp = freqz(sos,f,fs)*g;	% frequency response
m2_rsp = (f_rsp.*conj(f_rsp)).^2;	% squared for filtfilt
pow_psd = @(x) trapz(f(2:end),m2_rsp(2:end).*x(2:end));
pow_bg_filt = pow_psd(PSD_bg);
pow_sig_filt = pow_psd(PSD_smoo)-pow_bg_filt;
SNR_filt = pow_sig_filt/pow_bg_filt;    % SNR after filtering

%% Calculate detection threshold
auto_threshold = p.Results.auto_threshold;
threshold = p.Results.threshold;
if auto_threshold
    outliers = results.outliers;
    sig_range_i = results.sig_range_i;
    psd_cutoff = PSDcutoff(PSD_smoo,outliers,db2pow(results.tDB)*PSD_bg);
    pow_cutoff = trapz(f(sig_range_i(1):sig_range_i(2)), ...
        psd_cutoff(sig_range_i(1):sig_range_i(2)))-results.pow_bg;
    pow_cutoff_filt = pow_psd(psd_cutoff)-pow_bg_filt;
    prop = (pow_bg_filt+pow_cutoff_filt)/(pow_bg_filt+pow_sig_filt);
    LFP_amp_all = sort(LFP_amp_all);
    n = length(LFP_amp_all);
%     cumpow = cumsum(LFP_amp_all.^2);  % cumulative energy
    cumpow = cumsum(LFP_amp_all.^2)./(1:n)';  % cumulative power
    idx = max(find(cumpow>prop*cumpow(end),1)-1,1);
    threshold = LFP_amp_all(idx);
end
if isempty(threshold)
    ZS_threshold = p.Results.ZS_threshold;
    threshold = mean_amp+ZS_threshold*std_amp;
else
    threshold = max(threshold,0);
    ZS_threshold = (threshold-mean_amp)/std_amp;
end

%% Extract burst attributes
burst_seg = cell(nseg,1);
for i = 1:nseg
    burst_seg{i} = getDurations(LFP_amp{i},threshold,stop_perc);
    burst_seg{i} = getBurstAttr(burst_seg{i},LFP_seg{valid_seg_id(i)},fs,results.sig_range_f,Nfft);
    burst_seg{i}.pks_seg_id(:) = valid_seg_id(i);
    burst_seg{i}.bursts_seg_id(:) = valid_seg_id(i);
    burst_seg{i}.pkpha = angle(LFP_AS{i}(burst_seg{i}.tpks));
end
burst_seg = cell2mat(burst_seg);
stats = burst_seg(1);
if nseg>1
    FieldsSum = {'npeak','ndetect','nmerged','nvalid','runtime1','runtime2'};
    FieldsCat = {'peaks','tpks','pkpha','pks_seg_id','bursts_seg_id','duridx','AP','CN','BF'};
    for i = 1:length(FieldsSum)
        stats.(FieldsSum{i}) = sum([burst_seg.(FieldsSum{i})]);
    end
    for i = 1:length(FieldsCat)
        stats.(FieldsCat{i}) = vertcat(burst_seg.(FieldsCat{i}));
    end
end
AP_rate = stats.npeak/T_length;

if p.Results.fit_amp
    % Fit amplitude distribution
    nbins = 60;
    n = length(LFP_amp_all);
    maxamp = min(max(LFP_amp_all),10*mean_amp);
    LFP_amp_all = LFP_amp_all(LFP_amp_all<=maxamp);
    amp_edges = linspace(0,maxamp,nbins+1);
    amp_ctrs = bincenters(amp_edges);
    hist_amp = histcounts(LFP_amp_all,amp_edges);
    pdf_amp = hist_amp/n/binwidth(amp_edges);
    
    % Fit gamma distribution parameters for amplitude in signal component
    n = length(LFP_amp_all);
    if n>1e6
        rng(0); % random sample 1e6 points for performance
        LFP_amp_all = LFP_amp_all(randi(n,[1e6,1]));
    end
    gam_amp_tot = gamfit(LFP_amp_all);
    if p.Results.fit_amp_gamma
        tic;
        prec = 0.02;
        Grid = gamma_pdf_int_grid(gam_amp_tot(1),gam_amp_tot(2),prec,maxamp);
        warning('off','all');
        if p.Results.fit_amp_gamma == 1
            % Fit signal only
            mu_k = @(k) (pow_sig_filt*2/k/(k+1))^0.5;
            [gam_amp_sig,pci] = mle(LFP_amp_all,'pdf',@(x,k) ncxpdf_cond_gam(x,k,mu_k(k), ...
                pow_bg_filt^0.5,[],[],prec,Grid),'start',gam_amp_tot(1),'lowerbound',0.1);
            bg_pow_scale = 1;
            gam_amp_sig = [gam_amp_sig,mu_k(gam_amp_sig)];
        else
            % Fit background scale
            mu_k = @(k,s) ((pow_sig_filt+(1-s)*pow_bg_filt)*2/k/(k+1))^0.5;
            [fitparam,pci] = mle(LFP_amp_all,'pdf',@(x,k,s) ncxpdf_cond_gam(x,k,mu_k(k,s), ...
                (s*pow_bg_filt)^0.5,[],[],prec,Grid),'start',[gam_amp_tot(1),1], ...
                'lowerbound',[0.1,0.1],'upperbound',[Inf,1]);
            bg_pow_scale = fitparam(2);
            gam_amp_sig = [fitparam(1),mu_k(fitparam(1),bg_pow_scale)];
        end
        warning('on','all');
        runtime_fit = toc;
    end
end

%% Add results
varnames = {'unit','butter_order','auto_threshold','ZS_threshold','threshold', ...
    'stop_perc','freq_resolution','mean_amp','std_amp','AP_rate', ...
    'pow_sig_filt','pow_bg_filt','SNR_filt'};
if auto_threshold
    varnames = [varnames,'pow_cutoff','pow_cutoff_filt'];
end
if p.Results.fit_amp
    varnames = [varnames,'amp_edges','amp_ctrs','hist_amp','pdf_amp','gam_amp_tot'];
    if p.Results.fit_amp_gamma
        varnames = [varnames,'gam_amp_sig','bg_pow_scale','pci','runtime_fit'];
    end
end
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

%% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    if auto_threshold
        dispfcn('Auto threshold:');
    end
    dispfcn(num2str([threshold,ZS_threshold],['Detection threshold = %.5g ',unit,' (%.2f Z-score)']));
    dispfcn(num2str(AP_rate,'Total amplitude peak (AP) rate = %.2f Hz'));
    dispfcn(num2str(stats.ndetect/T_length,'Detected AP rate = %.2f Hz --- 100%%'));
    dispfcn(num2str([stats.nvalid/T_length,stats.nvalid/stats.ndetect*100],...
        'Valid detection rate = %.2f Hz --- %.2f%%'));
    deducted = stats.ndetect-stats.nmerged;
    dispfcn(num2str([deducted/T_length,deducted/stats.ndetect*100],...
        'AP rate deducted by merging = %.2f Hz --- %.2f%%'));
    deducted = stats.nmerged-stats.nvalid;
    dispfcn(num2str([deducted/T_length,deducted/stats.ndetect*100],...
        'AP rate deducted by failed frequency estimation = %.2f Hz --- %.2f%%'));
end

end