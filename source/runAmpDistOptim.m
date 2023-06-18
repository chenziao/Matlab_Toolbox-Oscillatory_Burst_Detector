function out = runAmpDistOptim(SynP,PSOparam,div_type,fit_AP,rerun,plotdist)
% SynP is an instance of "SynthParam" class.
% PSOparam is the structure containing PSO parameters.
% div_type is the amplitude distribution divergence type used for cost function.
% fit_AP is a flag to use amplitude peak as target for fitting. defualt: false, using amplitude.
% rerun is a flag to use background power scale as a parameter for optimization. default: false
% plotdist is the flag for plotting distribution after each evaluation. default: false
if nargin<4 || isempty(fit_AP),	fit_AP = false;	end
if nargin<5,	rerun = false;	end
if nargin<6,	plotdist = false;	end

% Determine parameter boundaries
amp_dist_type = SynP.amp_dist_type;
switch amp_dist_type
    case 'gamma'
        bd_scl = [1,15]/5;	% Upper/lower bound scaling factor
        gam_AP = SynP.rdat.gam_AP;
        LUB = [0.4,gam_AP(1)*bd_scl(2);gam_AP(2)*bd_scl];
        LUB(2,1) = LUB(2,1)/bd_scl(2);
        % LUB(2,2) = max(LUB(2,2),gam_AP(2)/LUB(1,1));
        LB = [1,LUB(2,1)];
    case 'lognormal'
        bd_scl = [1,15]/5;	% Upper/lower bound scaling factor
        logamp = SynP.rdat.logn_AP(1:2);
        LUB = [logamp(1)+max(log(bd_scl),[-inf,logamp(2)^2]); ...
            sqrt(max(logamp(2)^2-log(bd_scl([2,1])),[0.05^2,0])) ];
        LB = LUB(:,1);
    case 'exponential'
        bd_scl = [1,15]/5;   % Upper/lower bound scaling factor
        LUB = SynP.rdat.exp_AP*bd_scl;
        LB = LUB(1);
end
if rerun
    LUB = [LUB;0.1,1];
end

% Prepare for generation
gen_bursts_noAmp(SynP,LB,false);
SynP.setup_eval_amp_dist(div_type,fit_AP,[],[],[],[],plotdist);
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

if rerun
    SynP.amp_param.(amp_dist_type) = optimal_parameter(1:end-1);
    SynP.bg_pow_scale.(amp_dist_type) = optimal_parameter(end);
else
    SynP.amp_param.(amp_dist_type) = optimal_parameter;
end
% SynP.record_div_val(popCosts(imax));
SynP.record_div_val(sum(wts.*popCosts)/sum(wts));
eval_amp_dist(optimal_parameter,SynP,true);

end