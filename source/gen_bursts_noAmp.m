function SynP = gen_bursts_noAmp( SynP, LB, verbose )
% SynP is an instance of "SynthParam" class.
% LB are the lowerbound of the amplitude distribution parameters.
% verbose is the flag for displaying details such as progress. default: true
if ~isa(SynP,'SynthParam')
    error('First input must be an object of "SynthParam" class.');
end
if nargin<2
    LB = [];
elseif ~isempty(LB)
    switch SynP.amp_dist_type
        case 'exponential'
            if numel(LB)~=1,	error('Exponential distribution requires 1 parameter for the lowerbound.');	end
        case 'lognormal'
            if numel(LB)~=2,	error('Lognormal distribution requires 2 parameters for the lowerbound.');	end
        case 'gamma'
            if numel(LB)~=2,	error('Gamma distribution requires 2 parameters for the lowerbound.');	end
        otherwise
            error(['"',SynP.amp_dist_type,'" not available for amplitude distribution type.']);
    end
end
if nargin<3||~isscalar(verbose),	verbose = 1;	end
verbose = logical(verbose);

%% Estimate number of burst constrained by power
SynP.copula_pdf;
Y = SynP.cache.Xg(:,4)/SynP.dur_edge;
avg_wid = SynP.cache.PX'*Y*(SynP.width_edge/SynP.Dt);
if SynP.bl_power_match
    Y = SynP.bl_power_prop(exp(SynP.cache.Xg(:,3)),Y).*Y;
end
if isempty(LB)
    Nburst = inf;
else
    switch SynP.amp_dist_type
        case 'exponential'
            Y = expinv(normcdf(SynP.cache.Xg(:,1)),LB).^2.*Y;
        case 'lognormal'
            Y = exp((SynP.cache.Xg(:,1)*LB(2)+LB(1))*2).*Y;
        case 'gamma'
            Y = gaminv(normcdf(SynP.cache.Xg(:,1)),LB(1),LB(2)).^2.*Y;
    end
    Nburst = ceil(SynP.burst_energy/(SynP.cache.PX'*Y*SynP.unit_power));
end

try
    userview = memory;
    MemAvailable = userview.MemAvailableAllArrays;
    Windows = true;
catch
    MemAvailable = 4e9;
    Windows = false;
end
MemAllocate = max(MemAvailable-SynP.mem_presever,MemAvailable/2);
Nburst_memlim = floor(MemAllocate/(avg_wid+20)/8);
% Number of burst atoms to be generated
SynP.cache.Nburst = min([Nburst,Nburst_memlim,SynP.maxnumsampperdim^3]);
SynP.cache.Nburst = max(SynP.cache.Nburst,SynP.minnumburst);

if verbose
    if SynP.cache.Nburst==Nburst_memlim
        disp(num2str(MemAllocate/1e9,'Memory allocation limit %.2f GB is reached.'));
    elseif SynP.cache.Nburst==SynP.maxnumsampperdim^3
        disp('Maximum number of samples is reached.');
    elseif SynP.cache.Nburst==SynP.minnumburst
        disp('Minimum number of samples is reached.');
    end
    disp(num2str(SynP.cache.Nburst,'Number of burst atoms to be generated = %d.'));
end

%% Background energy
[PSD_bg,f,Nfft,bg_range_i] = SynP.psd_background;
if SynP.bl_power_match
    sig_range_f = SynP.rdat.sig_range_f;
    bg_range_i = [find(f>sig_range_f(1),1,'first'),find(f<sig_range_f(2),1,'last')];
    bg_range_i = [Nfft/2+2,bg_range_i(1):bg_range_i(2),Nfft/2+3];
    PSD_bg = [PSD_bg,SynP.rdat.sig_range_f];
    f = [f,sig_range_f];
end
SynP.cache.background_energy = trapz(f(bg_range_i),PSD_bg(bg_range_i))*SynP.Syn_len;

%% Generate Gaussian Copula and burst atom properties
% Set random seed for copula
if ~isempty(SynP.randseed),	rng(SynP.randseed+2);	end
if SynP.correlated
    SynP.cache.MVN = mvnrnd([0,0,0],SynP.sigcov,SynP.cache.Nburst);
else
    SynP.cache.MVN = randn([SynP.cache.Nburst,3]);
end
rndPhase = rand(SynP.cache.Nburst,1);

SynP.cache.burstCycNum = exp(SynP.cache.MVN(:,2)*SynP.sigma(2)+SynP.mu(2));
if SynP.empr_CN
    SynP.cache.burstCycNum = exp(interp1(SynP.rdat.cdf_log_CN,SynP.rdat.log_CNs,normcdf(SynP.cache.MVN(:,2))));
else
    SynP.cache.burstCycNum = exp(SynP.cache.MVN(:,2)*SynP.sigma(2)+SynP.mu(2));
end
if SynP.empr_BF
    SynP.cache.burstFreq = exp(interp1(SynP.rdat.cdf_log_BF,SynP.rdat.log_BFs,normcdf(SynP.cache.MVN(:,3))));
else
    SynP.cache.burstFreq = exp(SynP.cache.MVN(:,3)*SynP.sigma(3)+SynP.mu(3));
end
SynP.cache.Duration1sig = SynP.cache.burstCycNum./SynP.cache.burstFreq/SynP.dur_edge;
SynP.cache.wid = floor(SynP.width_edge/SynP.Dt*SynP.cache.Duration1sig)+1;
SynP.cache.ctr = ceil((SynP.cache.wid+1)/2);
sgl_ind = cumsum(SynP.cache.wid);
SynP.cache.sgl_ind = [[1;sgl_ind(1:end-1)+1],sgl_ind];
SynP.cache.sgl_burst = zeros(sum(SynP.cache.wid),1);

%% Generate unit amplitude atoms
npdf = 1000;
PDF = normpdf(linspace(-SynP.width_edge,SynP.width_edge,2*SynP.width_edge*npdf+1))'/normpdf(0);
if verbose
    disp(['Generating bursts ---',repmat(' ',1,19)]);
    bksp = repmat('\b',1,19);
    prog = -1;	prog_1000 = SynP.cache.Nburst/1000;	tic;
end
for i = 1:SynP.cache.Nburst
    tPts = (0:(SynP.cache.wid(i)-1))';
    SynP.cache.sgl_burst(SynP.cache.sgl_ind(i,1):SynP.cache.sgl_ind(i,2))...
        = PDF(round(SynP.Dt/SynP.cache.Duration1sig(i)*2*npdf*tPts+1))...
        .*sin(2*pi*(SynP.cache.burstFreq(i)*SynP.Dt*tPts+rndPhase(i)));
    
    if verbose
        progcurr = floor(i/prog_1000);
        if progcurr>prog
            prog = progcurr;
            fprintf([bksp,'%5.1f%% - %6.2f sec'],prog/10,toc);
        end
    end
end
if verbose,	fprintf([bksp,'%5.1f%% - %6.2f sec\n'],100,toc);	end

% Set random seed for bursts drawing
if ~isempty(SynP.randseed), rng(SynP.randseed+3);	end
SynP.cache.initrandstate = rng;	% record random state for inserting
SynP.cache.randstate = SynP.cache.initrandstate;

%% Check memory left
if verbose && Windows
    userview = memory;
    disp(num2str(SynP.mem_presever,'Memory preserved: %dB.'));
    disp(num2str(userview.MemAvailableAllArrays,'Memory available: %dB.'));
end

end