function [results,flag] = getPSDcomponents_fooof(results,fooof_results,freq_band,varargin)
p = inputParser;
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
validRangePair = @(x) validateattributes(x,{'numeric'},{'nonnegative','nondecreasing','numel',2});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',@(x) validateattributes(x,{'struct'},{'scalar'}));
addRequired(p,'fooof_results',@(x) validateattributes(x,{'struct'},{'scalar'}));
addRequired(p,'freq_band',validRangePair);  % signal frequency band
addOptional(p,'fit_range',[],validRangePair);   % frequency range used for fitting in FOOOF
addParameter(p,'tDB',1.0,validPositiveScalar);  % decibel threshold
% density of samples taken in log frequency for fitting (points per log Hz)
addParameter(p,'plot',false,validTF);  % plot PSD fit (set to 2 to plot autofit iterations)
addParameter(p,'display',false,validTF);
addParameter(p,'axes',[]);
addParameter(p,'msgfcn',[]);
parse(p,results,fooof_results,freq_band,varargin{:});

f = results.f;
aperiodic_params = fooof_results.aperiodic_params;
periodic_params = fooof_results.gaussian_params;
flag = -1;

fit_range = p.Results.fit_range;
freq_band(1) = max(freq_band(1),0);
freq_band(2) = min(freq_band(2),f(end));

%% fit curve
[PSD_smoo,PSD_fit,pe_fit] = getFOOOFcomponents(f,aperiodic_params,periodic_params);
fit_a = 10^aperiodic_params(1);
fit_b = aperiodic_params(end);

freq_range_i = [find(f>=freq_band(1),1),find(f<=freq_band(2),1,'last')];
if numel(freq_range_i)<2,	return; end

tDB = p.Results.tDB;
outliers = ind2seg(log(PSD_smoo)-log(PSD_fit)>log(db2pow(tDB)));
outliers = outliers(outliers(:,2)>=freq_range_i(1) & outliers(:,1)<=freq_range_i(2),:);

if isempty(fit_range)
    peak_params = fooof_results.peak_params;
    if size(peak_params,1)>0
        i_fit_range = [min(peak_params(:,1)-peak_params(:,3)),max(peak_params(:,1)+peak_params(:,3))];
        i_fit_range = [max(i_fit_range(1),f(2)),min(i_fit_range(2),f(end))];
        i_fit_range = [find(f<=i_fit_range(1),1,'last'),find(f>=i_fit_range(2),1,'first')];
    else
        i_fit_range = [2,numel(f)];
    end
else
    i_fit_range = [find(f(2:end)<=fit_range(1),1,'last')+1,find(f>=fit_range(2),1,'first')];
end

%% SNR
flag = 0;
if size(outliers,1)>0 && outliers(1)<outliers(end)
    flag = 1;
    sig_range_i = outliers([1,end]);
    sig_range_f = f(sig_range_i)';
    [~,troughs] = findpeaks(-pe_fit);
    intsec1 = troughs(find(troughs<sig_range_i(1),1,'last'));
    intsec2 = troughs(find(troughs>sig_range_i(2),1,'first'));
    if isempty(intsec1) || isempty(intsec2)
        psd_dif = results.PSD_smoo-PSD_fit;
        if isempty(intsec1)
            intsec1 = find(psd_dif(1:sig_range_i(1)-1)<0,1,'last')+1;
            if isempty(intsec1)
                [~,intsec1] = min(psd_dif(1:sig_range_i(1)-1));
                if isempty(intsec1),    intsec1 = 0;    end
            end
            intsec1 = max(intsec1,3);
        end
        if isempty(intsec2)
            intsec2 = find(psd_dif(sig_range_i(2)+1:end)<0,1,'first')-1;
            if isempty(intsec2)
                [~,intsec2] = min(psd_dif(sig_range_i(2)+1:end));
                if isempty(intsec2),    intsec2 = 0;    end
            end
            intsec2 = min(intsec2+sig_range_i(2),numel(f)-1);
        end
    end
    bg_range_i = intsec1:intsec2;
    bg_range_f = f([intsec1,intsec2])';
    i_fit_range = [min(i_fit_range(1),intsec1),max(i_fit_range(2),intsec2)];
    results.PSD_smoo = PSD_smoo;
    
    pow_psd = @(x) trapz(f(sig_range_i(1):sig_range_i(2)),x(sig_range_i(1):sig_range_i(2)));
    pow_bg = pow_psd(PSD_fit);
    pow_sig = pow_psd(PSD_smoo)-pow_bg;
    SNR = pow_sig/pow_bg;
else
    sig_range_i = [];    sig_range_f = [];
    bg_range_i = [];     bg_range_f = [];
    pow_bg = [];    pow_sig = [];   SNR = [];
end

%% Add results
varnames = {'tDB','freq_band','fit_range','i_fit_range','outliers', ...
    'sig_range_i','sig_range_f','bg_range_i','bg_range_f','fit_a','fit_b', ...
    'aperiodic_params','periodic_params','pow_bg','pow_sig','SNR'};
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

%% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    dispfcn(num2str(tDB,'Decibel threshold: %.2f dB'));
    if ~isempty(fit_range)
        dispfcn(['Frequency range for fitting: (',num2str(f(i_fit_range)','%.1f,%.1f'),')Hz']);
    end
    dispfcn(['Aperiodic parameters: ',num2str(aperiodic_params,'%.3f ')]);
    if flag>0
        dispfcn(['Signal frequency range: (',num2str(sig_range_f,'%.1f,%.1f'),')Hz']);
        if sig_range_f(1)<freq_band(1) || sig_range_f(2)>freq_band(2)
            dispfcn('Warning: The signal frequency range found exceeds the signal frequency band.');
        end
        dispfcn(['Background frequency range for substitution: (',num2str(bg_range_f,'%.1f,%.1f'),')Hz']);
        dispfcn(['Signal Noise Ratio: ',num2str(SNR,'%.2f '),' (',num2str(pow2db(SNR),'%.2f '),' dB)']);
    else
        dispfcn('No significant signal found within the signal frequency band.');
    end
end
if p.Results.plot
    PSDfit_plot(results,flag,p.Results.axes);
end

end