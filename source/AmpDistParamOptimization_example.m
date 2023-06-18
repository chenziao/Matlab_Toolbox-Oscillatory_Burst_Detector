close all;
clearvars -except working_directory;

%% Load
if ~exist('working_directory','var')
    working_directory = pwd;
end
[filename,working_directory] = uigetfile('*.mat','Select Input Data',working_directory);
if isequal(filename,0), clear working_directory;	return;	end
[~,fname] = fileparts(filename);
source_dir = pwd;
cd(working_directory);
addpath(source_dir);

%% General parameters settings
Amp_Types = {'gamma','lognormal','exponential'};
nPops = [12,12,8];
VarDisps = {[1,2,3],[1,2,3],[1,0,2]};

%% Check selected file
div_cutoff = 0.04;  % cutoff of KL/JS divergence to determine whether rerun
min_num_sample = 1e5;   % minimum number of samples required for rerunning optimization
data = load(filename);
if isfield(data,'SynParam')
    SynParam = data.SynParam;
	 %% Optimization settings
    [dist_type,avail_types] = SynParam.best_amp_dist_type;
    fit_AP = SynParam.eval_param.fit_AP;
    selected = cellfun(@(x) find(contains(Amp_Types,x)),avail_types);
    Div_Type = SynParam.eval_param.div_type;
    if SynParam.div_value.(dist_type)<=div_cutoff
        disp('Distribution matched well after first optimization.');	return;
    end
    disp('Distribution did not match well after first optimization. ');
    SynParam.loadrealdata;
    if fit_AP
        num_sample = round(SynParam.rdat.T_length*SynParam.rdat.AP_rate);
    else
        num_sample = round(SynParam.rdat.T_length*SynParam.rdat.fs);
    end
    if num_sample<min_num_sample
        disp('Number of samples is too few. Pass the second optimization.');	return;
    end
    disp('Running the second optimization which will scale down background power.');
    rerun = true;
else
    rerun = false;
    %% Optimization settings
    fit_AP = 1;   % whether fit amplitude peak (if not, fit amplitude)
    selected = 1:3;
    % selected = 1;
    Div_Type = 'KL';
    %% Setup
    SynParam = SynthParam;
    SynParam.filepath = filename;
    SynParam.randseed = 0;
    SynParam.mem_presever = 1e9;
    SynParam.empr_CN = true;
    SynParam.empr_BF = true;
    SynParam.syn_len = 1000;
    SynParam.loadrealdata;
    SynParam.fs = SynParam.rdat.fs;
end
SynParam.setup(1);
cd(source_dir);
clearvars data;

%% Generate noise trace
SynParam.gen_background;

%% PSO algorithm parameters
params.nPop = 8;            % Population Size (Swarm Size)
params.MaxIt = 30;          % Maximum Number of Iterations
params.MinIt = 15;          % Minimum Number of Iterations
params.MaxStall = 5;        % Terminate if change in best cost is ...
params.FunTol = 1e-5;       % less than FunTol over MaxStall iterations
params.withNeighbor = true; % Flag for using Neiborhood Operator
params.NeighborSize = [.3,1];	% Neighborhood Size (<=1 fraction, >1 number), [initial, final]
params.NeighborPeriod = 3;	% Number of iterations last for the same group of neighborhood
params.w = 1;               % Intertia Coefficient
params.wdamp = 0.92;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.InitVel = true;      % Random Non-Zero Initial Velocity
params.randseed = SynParam.randseed+4;	% Random seed
% PSO analysis parameters
params.ShowIterInfo = true;	% Flag for Showing Iteration Informatin
params.Visualize = true;

%% Run
for i = selected
    SynParam.amp_dist_type = Amp_Types{i};
    disp(['Optimizing with ',Amp_Types{i},' distribution:']);
    % PSO analysis parameters
    if rerun
        params.nPop = round(nPops(i)*1.5);
        params.VarDisp = VarDisps{i}([1,3]);
    else
        params.nPop = nPops(i);
        params.VarDisp = VarDisps{i}([1,2]);
    end
    tic;
    out = runAmpDistOptim(SynParam,params,Div_Type,fit_AP,rerun,0);
    toc
    % Results
    disp(num2str(out.FinalIteration,'Done after %d iterations.'));
    disp(num2str(out.exitflag,'Exitflag: %d'));
    switch out.exitflag
        case 1, disp('Maximum number of stall iterations reached.');
        case -1,    disp('Maximum number of iterations reached.');
    end
    if rerun
        disp(num2str(SynParam.bg_pow_scale.(Amp_Types{i}),'Background power scale: %.3f '));
    end
    disp(num2str(SynParam.div_value.(Amp_Types{i}),['Population ',Div_Type,' divergence = %f']));
    disp(num2str(SynParam.ks_test_pvalue.(Amp_Types{i}),'Kolmogorov-Smirnov test p-value = %f'));
    disp(['Optimized parameters: ',num2str(SynParam.amp_param.(Amp_Types{i}),'%.2f ')]);
    disp(' ');
    % Plot distribution
    SynParam.plot_eval_dist(Amp_Types{i});
end
if SynParam.div_value.(SynParam.best_amp_dist_type)>div_cutoff
    disp('Distribution did not match well.');
    if ~rerun
        disp('Rerun the optimization.');
    end
    disp(' ');
end

%% Save
if rerun
    fname = strrep(fname,'_SynParam','_SynParam_Rerun');
else
    fname = strrep(fname,'_Charac','_SynParam');
end
[filename,path] = uiputfile('*.mat','Save Results',fullfile(working_directory,fname));
if ~isequal(filename,0)
    working_directory = path;
    save(fullfile(working_directory,filename),'SynParam');
end
