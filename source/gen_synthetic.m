function [SynData,BurstAtoms] = gen_synthetic( SynP, DistP, verbose )
% Generate synthetic data returned in a structure "SynData".
% SynP is an instance of "SynthParam" class.
% distP is a vector that contains the parameters of the distribution.
% default: optimized parameters recorded in "SynP".
% verbose is the flag for displaying details such as progress. default: true
if ~isa(SynP,'SynthParam')
    error('First input must be an object of "SynthParam" class.');
end
if (nargin<2 || isempty(DistP)) && isfield(SynP.amp_param,SynP.amp_dist_type)
    DistP = SynP.amp_param.(SynP.amp_dist_type);
end
switch SynP.amp_dist_type
    case 'exponential'
        if numel(DistP)~=1, error('Exponential distribution requires 1 parameter.');	end
    case 'lognormal'
        if numel(DistP)~=2, error('Lognormal distribution requires 2 parameters.');	end
    case 'gamma'
        if numel(DistP)~=2,	error('Gamma distribution requires 2 parameters.');	end
    otherwise
        error(['"',SynP.amp_dist_type,'" not available for amplitude distribution type.']);
end
if nargin<3||~isscalar(verbose),	verbose = 1;	end
verbose = logical(verbose);

%% Generate noise trace
SynData.BackgroundTrace = SynP.gen_background;

%% Estimate number of burst constrained by power
[Xg,PX] = SynP.copula_pdf;
Y = Xg(:,4)*SynP.width_edge/SynP.dur_edge;
avg_wid = PX'*Y/SynP.Dt;
switch SynP.amp_dist_type
    case 'exponential'
        Y = expinv(normcdf(Xg(:,1)),DistP).^2.*Y;
    case 'lognormal'
        Y = exp((Xg(:,1)*DistP(2)+DistP(1))*2).*Y;
    case 'gamma'
        Y = gaminv(normcdf(Xg(:,1)),DistP(1),DistP(2)).^2.*Y;
end
nburst = round(SynP.burst_energy/(PX'*Y*SynP.unit_power));

try
    userview = memory;
    MemAvailable = userview.MemAvailableAllArrays;
catch
    MemAvailable = 4e9;
end
MemAllocate = max(MemAvailable-SynP.mem_presever,MemAvailable/2);
Nburst_memlim = floor(MemAllocate/(avg_wid+20)/8);
% Number of burst atoms to be generated
Nburst = min([nburst,Nburst_memlim,SynP.maxnumsampperdim^3]);

if verbose,
    if Nburst==Nburst_memlim
        disp(num2str(MemAllocate/1e9,'Memory allocation limit %.2f GB is reached.'));
    elseif Nburst==SynP.maxnumsampperdim^3
        disp('Maximum number of samples is reached.');
    end
    disp(num2str(Nburst,'Number of burst atoms to be generated = %d.'));
    disp(num2str(nburst/SynP.Syn_len,'Rate of burst atoms = %.2f Hz.'));
end

%% Generate Gaussian Copula and burst atom properties
% Set random seed for copula
if ~isempty(SynP.randseed),	rng(SynP.randseed+2);	end
if SynP.correlated
    MVN = mvnrnd([0,0,0],SynP.sigcov,Nburst);
else
    MVN = randn([Nburst,3]);
end
rndPhase = rand(Nburst,1);

switch SynP.amp_dist_type
    case 'exponential'
        burstAmp = expinv(normcdf(MVN(:,1)),DistP);
    case 'lognormal'
        burstAmp = exp(MVN(:,1)*DistP(2)+DistP(1));
    case 'gamma'
        burstAmp = gaminv(normcdf(MVN(:,1)),DistP(1),DistP(2));
end
if SynP.empr_CN
    burstCycNum = exp(interp1(SynP.rdat.cdf_cyc,SynP.rdat.log_cyc,normcdf(MVN(:,2))));
else
    burstCycNum = exp(MVN(:,2)*SynP.sigma(2)+SynP.mu(2));
end
if SynP.empr_BF
    burstFreq = exp(interp1(SynP.rdat.cdf_frq,SynP.rdat.log_frq,normcdf(MVN(:,3))));
else
    burstFreq = exp(MVN(:,3)*SynP.sigma(3)+SynP.mu(3));
end
Duration1sig = burstCycNum./burstFreq/SynP.dur_edge;

if SynP.power_match
    cumpow = cumsum(SynP.unit_power*SynP.width_edge*Duration1sig.*burstAmp.^2);
    nb = find(cumpow>=mod(SynP.burst_energy,cumpow(end)),1,'first');
    M = floor(SynP.burst_energy/cumpow(end));
    nburst = M*Nburst+nb;
else
    nb = mod(nburst-1,Nburst)+1;
    M = ceil(nburst/Nburst)-1;
end
Ibur = [repmat(1:Nburst,[1,M]),1:nb];
SynData.natom = nburst;

%% Generate and insert burst atoms
% Set random seed for inserting
if ~isempty(SynP.randseed), rng(SynP.randseed); end

npdf = 1000;
PDF = normpdf(linspace(-SynP.width_edge,SynP.width_edge,2*SynP.width_edge*npdf+1))'/normpdf(0);

SynData.SignalTrace = zeros(SynP.NT,1);
wid = floor(SynP.width_edge/SynP.Dt*Duration1sig(Ibur))+1;
burstSeq = sort(randi(SynP.NT,[nburst,1]));
t_insert = burstSeq-ceil((wid+1)/2);	% time index right before inserting
trueind = [t_insert+1,t_insert+wid];
ind_cmpr = [max(-t_insert,0),min(SynP.NT-trueind(:,2),0)];
trueind = trueind+ind_cmpr;
% durind = Duration1sig/2/SynP.Dt;
% durind = [ceil(durind*(SynP.width_edge-SynP.dur_edge))+t_insert+1,...
%     floor(durind*(SynP.width_edge+SynP.dur_edge))+t_insert+1];

if verbose
    disp(['Generating bursts ---',repmat(' ',1,19)]);
    bksp = repmat('\b',1,19);
    prog = -1;	prog_1000 = nburst/1000;	tic;
end
for i = 1:nburst
    I = Ibur(i);
    tPts = (ind_cmpr(i,1):wid(i)-1+ind_cmpr(i,2))';
    envelope = burstAmp(I)*PDF(round(SynP.Dt/Duration1sig(I)*2*npdf*tPts+1));
    sgl_burst = envelope.*sin(2*pi*(burstFreq(I)*SynP.Dt*tPts+rndPhase(I)));
    indList = trueind(i,1):trueind(i,2);
    SynData.SignalTrace(indList) = SynData.SignalTrace(indList)+sgl_burst;
    
    if verbose
        progcurr = floor(i/prog_1000);
        if progcurr>prog
            prog = progcurr;
            fprintf([bksp,'%5.1f%% - %6.2f sec'],prog/10,toc);
        end
    end
end
if verbose,	fprintf([bksp,'%5.1f%% - %6.2f sec\n'],100,toc);	end

% Collect Burst Atom Info
if nargout>1
    BurstAtoms = struct('AP',burstAmp,'CN',burstCycNum,'BF',burstFreq,'Ind',trueind);
end

end
