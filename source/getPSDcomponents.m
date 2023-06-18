function [results,flag] = getPSDcomponents(results,freq_band,varargin)
p = inputParser;
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
validRangePair = @(x) validateattributes(x,{'numeric'},{'nonnegative','nondecreasing','numel',2});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',@(x) validateattributes(x,{'struct'},{'scalar'}));
addRequired(p,'freq_band',validRangePair);  % signal frequency band
addOptional(p,'fit_range',[],validRangePair);   % frequency range for fitting
addParameter(p,'tDB',1.0,validPositiveScalar);  % decibel threshold
% density of sample frequency points in natural logarithm, neper (Np)
addParameter(p,'sample_density',50,validPositiveScalar);
addParameter(p,'autofit',false,validTF);    % find fit range automatically (set to 2 to enable tightmode)
addParameter(p,'plot',false,validTF);  % plot PSD fit (set to 2 to plot autofit iterations)
addParameter(p,'display',false,validTF);
addParameter(p,'axes',[]);  
addParameter(p,'msgfcn',[]);
parse(p,results,freq_band,varargin{:});

f = results.f;
PSD_smoo = results.PSD_smoo;
flag = -1;

fit_range = p.Results.fit_range;
autofit = p.Results.autofit || isempty(fit_range);
freq_band(1) = max(freq_band(1),0);
freq_band(2) = min(freq_band(2),f(end));

%% get PSD samples
if autofit
    fit_range = [f(2),f(end)];
else
    fit_range = max(fit_range,f(2));
end
i_fit = f>=fit_range(1) & f<=fit_range(2);
if isempty(i_fit),  return;	end
f_log = log(f(i_fit));	% find evenly distributed samples in log scale
n_fit = ceil((f_log(end)-f_log(1))*p.Results.sample_density);
f_pts = exp(linspace(f_log(1),f_log(end),n_fit));
i_fit = unique(discretize(f_pts,[0;f+binwidth(f)/2]));
if numel(i_fit)<3,  return;	end

%% fit curve
freq_range_i = [find(f(i_fit)>=freq_band(1),1),find(f(i_fit)<=freq_band(2),1,'last')];
if numel(freq_range_i)<2,	return; end
tDB = p.Results.tDB;
plot_iter = p.Results.plot>1;
if autofit
    [fit_a,fit_b,outliers,i_fit_range] = PSDfit(f(i_fit),PSD_smoo(i_fit), ...
        freq_range_i,tDB,'autofit',1,'tightmode',p.Results.autofit>1,'plot',plot_iter);
    i_fit_range = i_fit(i_fit_range);
else
    [fit_a,fit_b,outliers] = PSDfit(f(i_fit),PSD_smoo(i_fit), ...
        freq_range_i,tDB,'plot',plot_iter);
    i_fit_range = i_fit([1,end]);
end
PSD_fit = fit_a*f.^fit_b;
outliers = i_fit(outliers);
i_fit = i_fit_range(1):i_fit_range(2);
outliers_fine = i_fit(ind2seg(log(PSD_smoo(i_fit))-log(PSD_fit(i_fit))>log(db2pow(tDB))));
k = 1;
for i = 1:size(outliers,1)
    for j = k:size(outliers_fine,1)
        if outliers(i,1)>outliers_fine(j,1) && outliers(i,1)<=outliers_fine(j,2)
            outliers(i,1) = outliers_fine(j,1);
        end
        if outliers(i,2)>=outliers_fine(j,1) && outliers(i,2)<outliers_fine(j,2)
            outliers(i,2) = outliers_fine(j,2);
            k = j+1;    break;
        end
    end
end

%% SNR
flag = 0;
if size(outliers,1)>0 && outliers(1)<outliers(end)
    flag = 1;
    sig_range_i = outliers([1,end]);
    sig_range_f = f(sig_range_i)';
    psd_dif = PSD_smoo-PSD_fit;
    intsec1 = find(psd_dif(i_fit_range(1):sig_range_i(1)-1)<0,1,'last')+1;
    intsec2 = find(psd_dif(sig_range_i(2)+1:i_fit_range(2))<0,1,'first')-1;
    if isempty(intsec1)
        [~,intsec1] = min(psd_dif(i_fit_range(1):sig_range_i(1)-1));
    end
    if isempty(intsec2)
        [~,intsec2] = min(psd_dif(sig_range_i(2)+1:i_fit_range(2)));
    end
    if isempty(intsec1)
        intsec1 = i_fit_range(1);   flag = -3;
    else
        intsec1 = intsec1+i_fit_range(1)-1;
    end
    if isempty(intsec2)
        intsec2 = i_fit_range(2);   flag = -3;
    else
        intsec2 = intsec2+sig_range_i(2);
    end
    bg_range_i = intsec1:intsec2;
    bg_range_f = f([intsec1,intsec2])';
    
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
    'pow_bg','pow_sig','SNR'};
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

%% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    dispfcn(num2str(tDB,'Decibel threshold: %.2f dB'));
    dispfcn(['Frequency range for fitting: (',num2str(f(i_fit_range)','%.1f,%.1f'),')Hz']);
    dispfcn(num2str([fit_a,-fit_b],'Fit curve: S_B(f) = %.3f/f^%.3f'));
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