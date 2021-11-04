function [results,LFP_seg,flag] = getPSDcomponents(LFP_seg,fs,freq_range,varargin)
p = inputParser;
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
validRangePair = @(x) validateattributes(x,{'numeric'},{'nonnegative','nondecreasing','numel',2});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'LFP_seg',@(x) validateattributes(x,{'cell','numeric'},{'vector'}));
addRequired(p,'fs',validPositiveScalar);
addRequired(p,'freq_range',validRangePair);
addOptional(p,'fit_range',[]);
% scaled decibel threshold, multiple of 1 sigma
addParamValue(p,'nDB',1.645,validPositiveScalar);
% window for pwelch (>= 8 seconds)
addParamValue(p,'nfft',[],@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
% smoothing window width (Hz)
addParamValue(p,'smooth_win',2,validPositiveScalar);
% density of samples taken in log frequency for fitting (points per log Hz)
addParamValue(p,'fit_sample_density',50,validPositiveScalar);
addParamValue(p,'autofit',false,validTF);	% find fit range automatically
addParamValue(p,'plot',false,validTF);  % plot PSD fit (set to 2 to plot autofit iterations)
addParamValue(p,'display',false,validTF);
addParamValue(p,'axes',[]);
addParamValue(p,'msgfcn',[]);
parse(p,LFP_seg,fs,freq_range,varargin{:});

fit_range = p.Results.fit_range;
if ~isempty(fit_range) && ( ~isnumeric(fit_range) || numel(fit_range)~=2 )
    fit_range = [];
end
autofit = p.Results.autofit || isempty(fit_range);
if isempty(p.Results.nfft)
    nfft = 2^nextpow2(8*fs);
else
    nfft = p.Results.nfft;
end
noverlap = floor(nfft/2);   % 50% overlap
[~,f] = pwelch(zeros(nfft,1),nfft,noverlap,nfft,fs);
freq_range(1) = max(freq_range(1),0);
freq_range(2) = min(freq_range(2),f(end));

% Preprocess parameters
flag = -2;
results = [];
if ~iscell(LFP_seg)
    LFP_seg = {LFP_seg};
end
for i = 1:numel(LFP_seg)
    if ~isvector(LFP_seg{i}) || ~isnumeric(LFP_seg{i})
        error('Expect input data to be a cell vector, each element in which is a numeric vector.');
    end
    if ~isreal(LFP_seg{i})
        error('Expect input data to be real valued.');
    end
    LFP_seg{i} = double(LFP_seg{i}(:));
end
LFP_seg = LFP_seg(:);
LFP_mean = sum(cellfun(@sum,LFP_seg))/sum(cellfun(@length,LFP_seg));
LFP_seg = cellfun(@(x) x-LFP_mean,LFP_seg,'UniformOutput',0);

% remove short segments
size_seg = reshape(cellfun(@length,LFP_seg),[],1);
validseg = size_seg>=nfft;
T_length = sum(size_seg(validseg))/fs;
valid_seg_id = find(validseg);
if isempty(valid_seg_id),	return;	end
flag = -1;

% calculate PSD
nwindows = floor((size_seg(validseg)-noverlap)/(nfft-noverlap));	% # of windows
pxxs = cell2mat(reshape(cellfun(@(x) pwelch(x,nfft,noverlap,nfft,fs),LFP_seg(validseg),'UniformOutput',0),1,[]));
PSD = pxxs*nwindows/sum(nwindows);
PSD(1) = 0;

% smooth PSD
span = find(f<p.Results.smooth_win,1,'last');
PSD_smoo = moving(PSD,span);

% get PSD samples
if autofit
    [~,f_i] = max(PSD_smoo);
    fit_range = [min(f(f_i),2),f(end)];
else
    fit_range = max(fit_range,f(2));
end
i_fit = f>=fit_range(1) & f<=fit_range(2);
if isempty(i_fit),  return;	end
f_log = log(f(i_fit));	% find evenly distributed samples in log scale
n_fit = ceil((f_log(end)-f_log(1))*p.Results.fit_sample_density);
f_pts = exp(linspace(f_log(1),f_log(end),n_fit));
[~,i_fit] = min(abs(bsxfun(@minus,f_pts,f)),[],1);
% i_fit = discretize(f_pts,[f(1);conv(f,[1,1]/2,'valid');f(end)]);    % better in newer matlab version
i_fit = unique(i_fit);
if numel(i_fit)<3,  return;	end

% fit curve
freq_range_i = [find(f(i_fit)>=freq_range(1),1),find(f(i_fit)<=freq_range(2),1,'last')];
if numel(freq_range_i)<2,	return; end
nDB = p.Results.nDB;
plot_iter = p.Results.plot>1;
if autofit
    [fit_a,fit_b,outliers,i_fit_range] = PSDfit(f(i_fit),PSD_smoo(i_fit), ...
        freq_range_i,nDB,'autofit',1,'plot',plot_iter);
    i_fit_range = i_fit(i_fit_range);
else
    [fit_a,fit_b,outliers] = PSDfit(f(i_fit),PSD_smoo(i_fit), ...
        freq_range_i,nDB,'plot',plot_iter);
    i_fit_range = i_fit([1,end]);
end
outliers = i_fit(outliers);
PSD_fit = fit_a*f.^fit_b;

% SNR
flag = 0;
if size(outliers,1)>0 && outliers(1)<outliers(end)
    flag = 1;
    bump_i = outliers([1,end]);
    bump_f = f(bump_i)';
    psd_dif = PSD_smoo-PSD_fit;
    intsec1 = find(psd_dif(i_fit_range(1):bump_i(1)-1)<0,1,'last');
    intsec2 = find(psd_dif(bump_i(2)+1:i_fit_range(2))<0,1,'first');
    if isempty(intsec1)
        [~,intsec1] = min(psd_dif(i_fit_range(1):bump_i(1)-1));
    end
    if isempty(intsec2)
        [~,intsec2] = min(psd_dif(bump_i(2)+1:i_fit_range(2)));
    end
    if isempty(intsec1)
        intsec1 = i_fit_range(1);   flag = -3;
    else
        intsec1 = intsec1+i_fit_range(1);
    end
    if isempty(intsec2)
        intsec2 = i_fit_range(2);   flag = -3;
    else
        intsec2 = intsec2+bump_i(2)-1;
    end
    fit_i = intsec1:intsec2;
    fit_f = f([intsec1,intsec2])';
    
    pow_psd = @(x) trapz(f(bump_i(1):bump_i(2)),x(bump_i(1):bump_i(2)));
    pow_bg = pow_psd(PSD_fit);
    pow_sig = pow_psd(PSD_smoo)-pow_bg;
    SNR = pow_sig/pow_bg;
else
    bump_i = [];    bump_f = [];
    fit_i = [];     fit_f = [];
    pow_bg = [];    pow_sig = [];   SNR = [];
    if autofit && f(i_fit_range(1))>=freq_range(1) && f(i_fit_range(2))<=freq_range(2)
        % new flag: Suggestion: Increase the outlier threshold.
    end
end

% Results
varnames = {'T_length','fs','valid_seg_id','f','PSD','PSD_smoo','nDB','freq_range','i_fit_range',...
    'outliers','bump_i','bump_f','fit_i','fit_f','fit_a','fit_b','pow_bg','pow_sig','SNR'};
results = eval(['struct(',strjoin(cellfun(@(x) ['''',x,''',',x],varnames,'UniformOutput',0),','),')']);

% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    if length(validseg)>1
        dispfcn(num2str([length(valid_seg_id),length(validseg)],...
            '%d out of %d segments are long enough to use.'));
    end
    dispfcn(num2str(nDB,'Outlier threshold: %.2f dB'));
    dispfcn(['Frequency range for fitting: (',num2str(f(i_fit_range)','%.1f,%.1f'),')Hz']);
    dispfcn(num2str([fit_a,-fit_b],'Fit curve: S(f) = %.3f/f^%.3f'));
    if flag>0
        dispfcn(['Frequency range of signal of interest: (',num2str(bump_f,'%.1f,%.1f'),')Hz']);
        if bump_f(1)<freq_range(1) || bump_f(2)>freq_range(2)
            dispfcn('Warning: The frequency range of the signal found exceeds the bound of interest.');
        end
        dispfcn(['Intersection points: (',num2str(fit_f,'%.1f,%.1f'),')Hz']);
        dispfcn(['Signal Noise Ratio: ',num2str(SNR,'%.2f '),' (',num2str(pow2db(SNR),'%.2f '),' dB)']);
    else
        dispfcn('No significant signal found within the frequency range of interest.');
    end
end
if p.Results.plot
    PSDfit_plot(results,flag,p.Results.axes);
end

end

function c = moving(y,span)
% moving average of the data.
y = y(:);
span = floor(span);
n = length(y);
span = min(span,n);
width = span-1+mod(span,2); % force it to be odd
if width==1, c = y;	return;	end
c = filter(ones(width,1)/width,1,y);
cbegin = cumsum(y(1:width-2));
cbegin = cbegin(1:2:end)./(1:2:(width-2))';
cend = cumsum(y(n:-1:n-width+3));
cend = cend(end:-2:1)./(width-2:-2:1)';
c = [cbegin;c(width:end);cend];
end