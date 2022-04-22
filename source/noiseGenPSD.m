function y = noiseGenPSD(N,PSD,f,filt_n,source_type,filter_check,pow_match_freq)
% Use fir2 filter to generate noise given PSD
% N - number of samples to be generated
% PSD - One-sided PSD (Hz^-1) array. Containing total power in one side.
% f - Corresponding frequency (Hz) array (must be in ascending order). default: 0-0.5 Hz evenly spaced.
% Sampling frequency fs is assumed to be 2 times the maximun frequency in f.
% filt_n - order of fir2 filter (must be even valued). default: 500.
% source_type - distribution type of noise to be filtered. 'Gaussian'(default),'Uniform'.
% filter_check - Plot frequency response and compare generated PSD with desired PSD if true. default: false.
% pow_match_freq - the frequency range where generated power is forced to match the desired. default: [2*fs/filt_n,f(end)]

nPSD = numel(PSD);
PSD = reshape(PSD,1,nPSD);
if nargin < 3 || isempty(f)
    f = linspace(0,0.5,nPSD);
else
    if numel(f)~=nPSD
        error('The length of PSD and f must agree.');
    end
    f = reshape(f,1,nPSD);
    if f(1)~=0
        f = [0,f];
        PSD = [0,PSD];
        nPSD = nPSD+1;
    end
end
fs = 2*f(end);
if nargin < 4 || isempty(filt_n)
    filt_n = numel(f)-1;
end
filt_n = max(ceil(filt_n/2)*2,6);
if nargin < 5 || ~ischar(source_type)
    source_type = 'Gaussian';
end
if nargin < 6 || isempty(filter_check)
    filter_check = false;
end
if nargin < 7 || isempty(pow_match_freq) || sum(f>=pow_match_freq(1)&f<=pow_match_freq(2))<2
    pow_match_freq = [0,f(end)]+fs/2/filt_n*[1,-1];
end

m = sqrt(PSD);
b = fir2(filt_n,f*2/fs,m);
f_rsp = freqz(b,1,f,fs);
m2_rsp = f_rsp.*conj(f_rsp);
ind = f>=pow_match_freq(1) & f<=pow_match_freq(2);
power = trapz(f(ind),PSD(ind))*trapz(f,m2_rsp)/trapz(f(ind),m2_rsp(ind));

if filter_check
    try
        figure;
        freqz(b);
    catch
        close(gcf);
    end
end

switch source_type
    case 'Uniform'
        x = rand(N+filt_n,1);
    case 'Laplace'
        x = 2*rand(N+filt_n,1)-1;
        x = -sign(x).*log(1-abs(x));
    otherwise	% Gaussian
        x = randn(N+filt_n,1);
end

y = filter(b,1,x);
y = sqrt(power)*zscore(y(filt_n+1:end));

if filter_check
    nfft = (nPSD-1)*2;
    [pxx,ff] = pwelch(y,nfft,0,nfft,fs);
    
    figure;
    plot(ff,pow2db(pxx),f,pow2db(PSD));
    legend({'Generated PSD','Desired PSD'});
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
end

end