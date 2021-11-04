close all;
clearvars -except working_directory;

%% Load
if ~exist('working_directory','var')
    working_directory = fullfile('examples','example1');
end
[filename,working_directory] = uigetfile('*.mat','Select Input Data',working_directory);
if isequal(filename,0), clear working_directory;	return;	end
[~,fname] = fileparts(filename);
source_dir = pwd;
cd(working_directory);
addpath(source_dir);

%% Optimization settings
fit_ASAP = 1;   % Use AS-AP (default: 0 AS amplitude)
Amp_Types = {'exponential','gamma','lognormal'};
selected = [2,1,3];
nPops = [8,12,12];
VarDisps = {[1,0],[1,2],[1,2]};
Div_Type = 'KL';

%% Setup
SynParam = SynthParam;
SynParam.filepath = filename;
SynParam.randseed = 0;
SynParam.mem_presever = 1e9;
% SynParam.butter_order = 2;
SynParam.empr_CN = true;
SynParam.empr_BF = true;
SynParam.syn_len = 1000;
SynParam.loadrealdata;
SynParam.fs = SynParam.rdat.fs;
SynParam.setup;

cd(source_dir);

%% Generate noise trace
SynParam.gen_background;

%% PSO algorithm parameters
params.nPop = 8;            % Population Size (Swarm Size)
params.MaxIt = 30;          % Maximum Number of Iterations
params.MinIt = 15;          % Minimum Number of Iterations
params.MaxStall = 5;        % Terminate if change in best cost is ...
params.FunTol = 1e-5;       % less than FunTol over MaxStall iterations
params.withNeighbor = true; % Flag for using Neiborhood Operator
params.NeighborSize = [.3,1];	% Neighborhood Size (<=1 fraction, >1 number)
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
    params.nPop = nPops(i);
    params.VarDisp = VarDisps{i};
    tic;
    runAmpDistOptim(SynParam,params,Div_Type,fit_ASAP,0);
    toc
    % Results
    disp([Amp_Types{i},' parameters: ',num2str(SynParam.amp_param.(Amp_Types{i}),'%.2f ')]);
    disp(num2str(SynParam.div_value.(Amp_Types{i}),['Population ',Div_Type,' divergence = %f']));
    disp(' ');
    % Plot distribution
    SynParam.plot_eval_dist(Amp_Types{i});
end

%% Save
fname = strrep(fname,'_Charac','_SynParam');
[filename,path] = uiputfile('*.mat','Save Results',fullfile(working_directory,fname));
if ~isequal(filename,0)
    working_directory = path;
    save(fullfile(working_directory,filename),'SynParam');
end
