function [Div,h1,h2] = eval_amp_dist( DistP, SynP, save_dist, varargin )
% Evaluate the divergence between the amplitude distribution of the 
% real data and the synthetic data.
% distP is a vector that contains the parameters of the distribution.
% SynP is an instance of "SynthParam" class.
% save_dist is a flag indicating whether to save the distribution for this evaluation.
if ~isa(SynP,'SynthParam')
    error('First input must be an object of "SynthParam" class.');
end
amp_dist_type = SynP.amp_dist_type;
switch amp_dist_type
    case 'exponential'
        if numel(DistP)~=1, error('Exponential distribution requires 1 parameter.');	end
    case 'lognormal'
        if numel(DistP)~=2, error('Lognormal distribution requires 2 parameters.');	end
    case 'gamma'
        if numel(DistP)~=2, error('Gamma distribution requires 2 parameters.');	end
    otherwise
        error(['"',amp_dist_type,'" not available for amplitude distribution type.']);
end
if isempty(fieldnames(SynP.eval_param))
    error('Parameters need to be setup using function "setup_eval_amp_dist" of "SynthParam" class.');
end
if nargin<3 || isempty(save_dist)
    save_dist = false;
end

divfun = @(x) SynP.eval_param.divfun(x(SynP.eval_param.threshold_idx:end));
fit_ASAP = SynP.eval_param.fit_ASAP;
randburst = SynP.eval_param.randburst;
verbose = SynP.eval_param.verbose;
plotdist = SynP.eval_param.plotdist;

%% Determine subset of burst atoms to be inserted and their order
if randburst,	rng(SynP.buffer.randstate); end	% restore random state
if SynP.power_match
    if randburst,	Iperm  = randperm(SynP.buffer.Nburst);
    else	Iperm = 1:SynP.buffer.Nburst;	end
else
    switch amp_dist_type
        case 'exponential'
            Y = expinv(normcdf(SynP.buffer.Xg(:,1)),DistP).^2.*SynP.buffer.Xg(:,4);
        case 'lognormal'
            Y = exp((SynP.buffer.Xg(:,1)*DistP(2)+DistP(1))*2).*SynP.buffer.Xg(:,4);
        case 'gamma'
            Y = gaminv(normcdf(SynP.buffer.Xg(:,1)),DistP(1),DistP(2)).^2.*SynP.buffer.Xg(:,4);
    end
    nburst = ceil(SynP.burst_energy/(SynP.buffer.PX'*Y*(SynP.width_edge/SynP.dur_edge*SynP.unit_power)));
    nb = mod(nburst-1,SynP.buffer.Nburst)+1;
    M = ceil(nburst/SynP.buffer.Nburst)-1;
    NB = min(nburst,SynP.buffer.Nburst);
    if randburst,	Iperm  = randperm(SynP.buffer.Nburst,NB);
    else	Iperm = 1:NB;	end
end
% generate amplitudes
switch amp_dist_type
    case 'exponential'
        burstAmp = expinv(normcdf(SynP.buffer.MVN(Iperm,1)),DistP);
    case 'lognormal'
        burstAmp = exp(SynP.buffer.MVN(Iperm,1)*DistP(2)+DistP(1));
    case 'gamma'
        burstAmp = gaminv(normcdf(SynP.buffer.MVN(Iperm,1)),DistP(1),DistP(2));
end
if SynP.power_match
    cumpow = cumsum(SynP.unit_power*SynP.width_edge*SynP.buffer.Duration1sig.*burstAmp.^2);
    nb = find(cumpow>=mod(SynP.burst_energy,cumpow(end)),1,'first');
    M = floor(SynP.burst_energy/cumpow(end));
    nburst = M*SynP.buffer.Nburst+nb;
end

II = [repmat(1:SynP.buffer.Nburst,[1,M]),1:nb];
Ibur = Iperm(II);
if randburst,	SynP.buffer.randstate = rng;	end	% record random state for next inserting

%% Insert busrts
if ~isempty(SynP.randseed),	rng(SynP.randseed);	end

signalTrace = zeros(SynP.NT,1);
t_insert = sort(randi(SynP.NT,[nburst,1])) - SynP.buffer.ctr(Ibur);
trueind = [t_insert+1,t_insert+SynP.buffer.wid(Ibur)];
ind_cmpr = [max(-t_insert,0),min(SynP.NT-trueind(:,2),0)];
trueind = trueind+ind_cmpr;
Sgl_Ind = SynP.buffer.sgl_ind(Ibur,:)+ind_cmpr;
if verbose
    disp(['Inserting bursts ---',repmat(' ',1,17)]);
    bksp = repmat('\b',1,17);
    prog = -1;  prog_100 = nburst/100;	tic;
    for i = 1:nburst
        signalTrace(trueind(i,1):trueind(i,2)) = signalTrace(trueind(i,1):trueind(i,2)) + ...
            burstAmp(II(i))*SynP.buffer.sgl_burst(Sgl_Ind(i,1):Sgl_Ind(i,2));
        progcurr = floor(i/prog_100);
        if progcurr>prog
            prog = progcurr;
            fprintf([bksp,'%3.0f%% - %6.2f sec'],prog,toc);
        end
    end
else
    for i = 1:nburst
        signalTrace(trueind(i,1):trueind(i,2)) = signalTrace(trueind(i,1):trueind(i,2)) + ...
            burstAmp(II(i))*SynP.buffer.sgl_burst(Sgl_Ind(i,1):Sgl_Ind(i,2));
    end
end

%% Get AS amplitude
% LFP_filt = filtfilt(SynP.bfilt,SynP.afilt,SynP.buffer.backgroundTrace+signalTrace);
LFP_filt = filtfilt(SynP.sos,SynP.g,SynP.buffer.backgroundTrace+signalTrace);	% using sos
AS_amp = abs(hilbert(LFP_filt));
if verbose, fprintf([bksp,'%3.0f%% - %6.2f sec\n'],100,toc);    end

%% Evaluate distribution
% fit AS peaks
if fit_ASAP || save_dist || plotdist
    peaks = Findpeaks(AS_amp);
    hist_pks = hist(peaks,SynP.rdat.pk_amp);
end

% fit AS amplitude
if ~fit_ASAP || save_dist
    if SynP.NT>1e6
        AS_amp = AS_amp(randi(SynP.NT,[1e6,1]));
    end
    hist_amp = hist(AS_amp,SynP.rdat.as_amp);
end

if save_dist
    SynP.dist_snapshot.(amp_dist_type) = {hist_amp/numel(AS_amp)/binwidth(SynP.rdat.as_amp),...
        hist_pks/SynP.syn_len/binwidth(SynP.rdat.pk_amp)};
else
    if fit_ASAP
        Div = divfun(hist_pks);
    else
        Div = divfun(hist_amp);
    end
    if plotdist,	[h1,h2] = SynP.plot_peak_dist(hist_pks,varargin{:});	end
end

end
