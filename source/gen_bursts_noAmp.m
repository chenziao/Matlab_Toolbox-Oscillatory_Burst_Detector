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
Y = SynP.buffer.Xg(:,4)*SynP.width_edge/SynP.dur_edge;
avg_wid = SynP.buffer.PX'*Y/SynP.Dt;
if isempty(LB)
    Nburst = inf;
else
    switch SynP.amp_dist_type
        case 'exponential'
            Y = expinv(normcdf(SynP.buffer.Xg(:,1)),LB).^2.*Y;
        case 'lognormal'
            Y = exp((SynP.buffer.Xg(:,1)*LB(2)+LB(1))*2).*Y;
        case 'gamma'
            Y = gaminv(normcdf(SynP.buffer.Xg(:,1)),LB(1),LB(2)).^2.*Y;
    end
    Nburst = round(SynP.burst_energy/(SynP.buffer.PX'*Y*SynP.unit_power));
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
SynP.buffer.Nburst = min([Nburst,Nburst_memlim,SynP.maxnumsampperdim^3]);

if verbose
    if SynP.buffer.Nburst==Nburst_memlim
        disp(num2str(MemAllocate/1e9,'Memory allocation limit %.2f GB is reached.'));
    elseif SynP.buffer.Nburst==SynP.maxnumsampperdim^3
        disp('Maximum number of samples is reached.');
    end
    disp(num2str(SynP.buffer.Nburst,'Number of burst atoms to be generated = %d.'));
end

%% Generate Gaussian Copula and burst atom properties
% Set random seed for copula
if ~isempty(SynP.randseed),	rng(SynP.randseed+2);	end
if SynP.correlated
    SynP.buffer.MVN = mvnrnd([0,0,0],SynP.sigcov,SynP.buffer.Nburst);
else
    SynP.buffer.MVN = randn([SynP.buffer.Nburst,3]);
end
rndPhase = rand(SynP.buffer.Nburst,1);

SynP.buffer.burstCycNum = exp(SynP.buffer.MVN(:,2)*SynP.sigma(2)+SynP.mu(2));
if SynP.empr_CN
    SynP.buffer.burstCycNum = exp(interp1(SynP.rdat.cdf_cyc,SynP.rdat.log_cyc,normcdf(SynP.buffer.MVN(:,2))));
else
    SynP.buffer.burstCycNum = exp(SynP.buffer.MVN(:,2)*SynP.sigma(2)+SynP.mu(2));
end
if SynP.empr_BF
    SynP.buffer.burstFreq = exp(interp1(SynP.rdat.cdf_frq,SynP.rdat.log_frq,normcdf(SynP.buffer.MVN(:,3))));
else
    SynP.buffer.burstFreq = exp(SynP.buffer.MVN(:,3)*SynP.sigma(3)+SynP.mu(3));
end
SynP.buffer.Duration1sig = SynP.buffer.burstCycNum./SynP.buffer.burstFreq/SynP.dur_edge;
SynP.buffer.wid = floor(SynP.width_edge/SynP.Dt*SynP.buffer.Duration1sig)+1;
SynP.buffer.ctr = ceil((SynP.buffer.wid+1)/2);
sgl_ind = cumsum(SynP.buffer.wid);
SynP.buffer.sgl_ind = [[1;sgl_ind(1:end-1)+1],sgl_ind];
SynP.buffer.sgl_burst = zeros(sum(SynP.buffer.wid),1);

%% Generate unit amplitude atoms
npdf = 1000;
PDF = normpdf(linspace(-SynP.width_edge,SynP.width_edge,2*SynP.width_edge*npdf+1))'/normpdf(0);
if verbose
    disp(['Generating bursts ---',repmat(' ',1,19)]);
    bksp = repmat('\b',1,19);
    prog = -1;	prog_1000 = SynP.buffer.Nburst/1000;	tic;
end
for i = 1:SynP.buffer.Nburst
    tPts = (0:(SynP.buffer.wid(i)-1))';
    SynP.buffer.sgl_burst(SynP.buffer.sgl_ind(i,1):SynP.buffer.sgl_ind(i,2))...
        = PDF(round(SynP.Dt/SynP.buffer.Duration1sig(i)*2*npdf*tPts+1))...
        .*sin(2*pi*(SynP.buffer.burstFreq(i)*SynP.Dt*tPts+rndPhase(i)));
    
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
SynP.buffer.initrandstate = rng;	% record random state for inserting
SynP.buffer.randstate = SynP.buffer.initrandstate;

%% Check memory left
if verbose && Windows
    userview = memory;
    disp(num2str(SynP.mem_presever,'Memory preserved: %dB.'));
    disp(num2str(userview.MemAvailableAllArrays,'Memory available: %dB.'));
end

end