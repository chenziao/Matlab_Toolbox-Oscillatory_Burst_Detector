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
if isempty(fieldnames(SynP.eval_param))
    error('Parameters need to be setup using function "setup_eval_amp_dist" of "SynthParam" class.');
end
if nargin<3 || isempty(save_dist)
    save_dist = false;
end

bg_pow_scale = 1;
switch amp_dist_type
    case 'exponential'
        if numel(DistP)>1,	bg_pow_scale = DistP(2);	end
    case 'lognormal'
        if numel(DistP)>2,	bg_pow_scale = DistP(3);	end
    case 'gamma'
        if numel(DistP)>2,	bg_pow_scale = DistP(3);	end
end

divfun = @(x) SynP.eval_param.divfun(x(SynP.eval_param.threshold_idx:end));
fit_AP = SynP.eval_param.fit_AP;
randburst = SynP.eval_param.randburst;
verbose = SynP.eval_param.verbose;
plotdist = SynP.eval_param.plotdist;
if isfield(SynP.eval_param,'kstest')
    ks = logical(SynP.eval_param.kstest);
else
    ks = false;
end

%% Determine subset of burst atoms to be inserted and their order
if randburst,	rng(SynP.cache.randstate); end	% restore random state
burst_energy = SynP.burst_energy+(1-bg_pow_scale)*SynP.cache.background_energy;
if SynP.power_match
    if randburst,	Iperm  = randperm(SynP.cache.Nburst);
    else	Iperm = 1:SynP.cache.Nburst;	end
else
    Y = SynP.cache.Xg(:,4)/SynP.dur_edge;
    if SynP.bl_power_match
        Y = SynP.bl_power_prop(exp(SynP.cache.Xg(:,3)),Y).*Y;
    end
    switch amp_dist_type
        case 'exponential'
            Y = expinv(normcdf(SynP.cache.Xg(:,1)),DistP(1)).^2.*Y;
        case 'lognormal'
            Y = exp((SynP.cache.Xg(:,1)*DistP(2)+DistP(1))*2).*Y;
        case 'gamma'
            Y = gaminv(normcdf(SynP.cache.Xg(:,1)),DistP(1),DistP(2)).^2.*Y;
    end
    nburst = ceil(burst_energy/(SynP.cache.PX'*Y*SynP.unit_power));
    nb = mod(nburst-1,SynP.cache.Nburst)+1;
    M = ceil(nburst/SynP.cache.Nburst)-1;
    NB = min(nburst,SynP.cache.Nburst);
    if randburst,	Iperm  = randperm(SynP.cache.Nburst,NB);
    else	Iperm = 1:NB;	end
end
% generate amplitudes
switch amp_dist_type
    case 'exponential'
        burstAmp = expinv(normcdf(SynP.cache.MVN(Iperm,1)),DistP(1));
    case 'lognormal'
        burstAmp = exp(SynP.cache.MVN(Iperm,1)*DistP(2)+DistP(1));
    case 'gamma'
        burstAmp = gaminv(normcdf(SynP.cache.MVN(Iperm,1)),DistP(1),DistP(2));
end
if SynP.power_match
    if SynP.bl_power_match
        cumpow = cumsum(SynP.unit_power*SynP.cache.Duration1sig.*burstAmp.^2 ...
            .*SynP.bl_power_prop(SynP.cache.burstFreq,SynP.cache.Duration1sig));
    else
        cumpow = cumsum(SynP.unit_power*SynP.cache.Duration1sig.*burstAmp.^2);
    end
    nb = find(cumpow>=mod(burst_energy,cumpow(end)),1,'first');
    M = floor(burst_energy/cumpow(end));
    nburst = M*SynP.cache.Nburst+nb;
end

II = [repmat(1:SynP.cache.Nburst,[1,M]),1:nb];
Ibur = Iperm(II);
if randburst,	SynP.cache.randstate = rng;	end	% record random state for next inserting

%% Insert busrts
if ~isempty(SynP.randseed),	rng(SynP.randseed);	end

signalTrace = zeros(SynP.NT,1);
t_insert = sort(randi(SynP.NT,[nburst,1])) - SynP.cache.ctr(Ibur);
trueind = [t_insert+1,t_insert+SynP.cache.wid(Ibur)];
ind_cmpr = [max(-t_insert,0),min(SynP.NT-trueind(:,2),0)];
trueind = trueind+ind_cmpr;
Sgl_Ind = SynP.cache.sgl_ind(Ibur,:)+ind_cmpr;
if verbose
    disp(['Inserting bursts ---',repmat(' ',1,17)]);
    bksp = repmat('\b',1,17);
    prog = -1;  prog_100 = nburst/100;	tic;
    for i = 1:nburst
        signalTrace(trueind(i,1):trueind(i,2)) = signalTrace(trueind(i,1):trueind(i,2)) + ...
            burstAmp(II(i))*SynP.cache.sgl_burst(Sgl_Ind(i,1):Sgl_Ind(i,2));
        progcurr = floor(i/prog_100);
        if progcurr>prog
            prog = progcurr;
            fprintf([bksp,'%3.0f%% - %6.2f sec'],prog,toc);
        end
    end
else
    for i = 1:nburst
        signalTrace(trueind(i,1):trueind(i,2)) = signalTrace(trueind(i,1):trueind(i,2)) + ...
            burstAmp(II(i))*SynP.cache.sgl_burst(Sgl_Ind(i,1):Sgl_Ind(i,2));
    end
end

%% Get amplitude
if bg_pow_scale<1
    LFP = bg_pow_scale^.5*SynP.cache.backgroundTrace+signalTrace;
else
    LFP = SynP.cache.backgroundTrace+signalTrace;
end
amp = abs(hilbert(filtfilt(SynP.sos,SynP.g,LFP)));   % using sos
if verbose, fprintf([bksp,'%3.0f%% - %6.2f sec\n'],100,toc);    end

%% Evaluate distribution
% fit amplitude peaks
if fit_AP || save_dist || plotdist
    AP = Findpeaks(amp);
    if ~ks || save_dist || plotdist
        hist_AP = histcounts(AP,SynP.rdat.AP_edges);
    end
end

% fit amplitude
if ~fit_AP || save_dist
    if SynP.NT>1e6
        amp = amp(randi(SynP.NT,[1e6,1]));
    end
    if ~ks || save_dist || plotdist
        hist_amp = histcounts(amp,SynP.rdat.amp_edges);
    end
end

if ks || save_dist
    if fit_AP
        CDF = [SynP.rdat.AP_edges;0,cumsum(SynP.rdat.hist_AP)]';
        m = max(AP);
    else
        CDF = [SynP.rdat.amp_edges;0,cumsum(SynP.rdat.hist_amp)]';
        m = max(amp);
    end
    CDF(:,2) = CDF(:,2)/CDF(end,2);
    if m>CDF(end,1)
        CDF = [CDF;m,1];
    end
    if fit_AP
        [~,p,ksstat] = kstest(AP,'CDF',CDF);
    else
        [~,p,ksstat] = kstest(amp,'CDF',CDF);
    end
end

if save_dist
    SynP.dist_snapshot.(amp_dist_type) = {hist_amp/numel(amp)/binwidth(SynP.rdat.amp_edges),...
        hist_AP/SynP.syn_len/binwidth(SynP.rdat.AP_edges)};
    SynP.ks_test_pvalue.(amp_dist_type) = p;
else
    if ~ks
        if fit_AP
            Div = divfun(hist_AP);
        else
            Div = divfun(hist_amp);
        end
    else
        Div = ksstat;   % -log(p)
    end
    if plotdist,	[h1,h2] = SynP.plot_AP_dist(hist_AP,varargin{:});	end
end

end
