function out = runAmpDistOptim(SynP,PSOparam,div_type,fit_ASAP,plotdist)
% SynP is an instance of "SynthParam" class.
% PSOparam is the structure containing PSO parameters.
% div_type is the amplitude distribution type used for generation.
% fit_ASAP is a flag to use AS-AP as objective for fitting. defualt: false, using AS amplitude.
% plotdist is the flag for plotting distribution after each evaluation. default: false
if nargin<4 || isempty(fit_ASAP),	fit_ASAP = false;	end
if nargin<5,	plotdist = false;	end

% Determine parameter boundaries
amp_dist_type = SynP.amp_dist_type;
switch amp_dist_type
    case 'exponential'
        bd_scl = [1,10]/5;   % Upper/lower bound scaling factor
        LUB = SynP.rdat.expAP*bd_scl;
        LB = LUB(1);
    case 'gamma'
        bd_scl = [1,10]/5;	% Upper/lower bound scaling factor
        gamamp = SynP.rdat.gamAP;
        LUB = [0.4,gamamp(1)*bd_scl(2);gamamp(2)*bd_scl];
        LUB(2,1) = LUB(2,1)/bd_scl(2);
        LUB(2,2) = max(LUB(2,2),gamamp(1)/LUB(1,1));
        LB = [1,LUB(2,1)];
    case 'lognormal'
        bd_scl = [1,10]/5;	% Upper/lower bound scaling factor
        logamp = SynP.rdat.logAP(1:2);
        LUB = [logamp(1)+max(log(bd_scl),[-inf,logamp(2)^2]); ...
            sqrt(max(logamp(2)^2-log(bd_scl([2,1])),[0.05^2,0])) ];
        LB = LUB(:,1);
end

% Prepare for generation
gen_bursts_noAmp(SynP,LB,false);
SynP.setup_eval_amp_dist(div_type,fit_ASAP,[],[],[],[],plotdist);
SynP.reset_randstate;

% Problem definition
problem.CostFunction = @(x) eval_amp_dist(x,SynP,false);	% Cost Function
problem.LUbounds = LUB;

% Run Optimization
out = PSO(problem, PSOparam);

% Find solution
pop = out.pop;
popSol = cell2mat(arrayfun(@(x) x.Best.Position,pop,'UniformOutput',false));
popCosts = arrayfun(@(x) x.Best.Cost,pop);
[~,rank] = sort(popCosts);
npop = length(pop);
wts = zeros(npop,1);
wts(rank) = cos(linspace(0,pi,npop))+1;
try
    f = ksdensity(popSol,popSol,'weights',wts);	% bivariate
catch
    f = ksdensity(popSol(:,1),popSol(:,1),'weights',wts);	% univariate
end
[~,imax] = max(f);
optimal_parameter = popSol(imax,:);

SynP.amp_param.(amp_dist_type) = optimal_parameter;
% SynP.record_div_val(popCosts(imax));
SynP.record_div_val(sum(wts.*popCosts)/sum(wts));
eval_amp_dist(optimal_parameter,SynP,true);

end