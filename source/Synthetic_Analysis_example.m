close all;
clearvars -except working_directory;

if ~exist('working_directory','var')
    working_directory = pwd;
end
[filename,working_directory] = uigetfile('*.mat','Select Input Data',working_directory);
if isequal(filename,0),	clear working_directory;	return;	end
filename = fullfile(working_directory,filename);
source_dir = pwd;
cd(working_directory);
addpath(source_dir);

%% Load synthetic data
load(filename);
SynParam.loadrealdata;
cd(source_dir);

if ~exist('SynData','var')
    SynParam.syn_len = 5000;
    select_type = 'gamma';
    if isempty(select_type)
        dist_type = SynParam.best_amp_dist_type;
        if isempty(dist_type)
            return;
        end
        SynParam.amp_dist_type = dist_type;
    else
        SynParam.amp_dist_type = select_type;
    end
    SynParam.setup;
    SynData = gen_synthetic(SynParam);
end
disp(['Amplitude distribution type being used: ',SynParam.amp_dist_type]);

%% Primary Process
SynStats = getSyntheticStats(SynParam,SynData);

%% Preselected bounds
BoundStats = getPreselectedBounds(SynParam,SynStats);
disp(num2str(BoundStats.nsigma_bound,[char(952),'_L = %.2f',char(963), ...
    ', ',char(952),'_U = %.2f',char(963)]));

%% Detection
detection_pk = getDetection(SynStats,BoundStats,'Detect_ASAP',true);
detection_amp = getDetection(SynStats,BoundStats,'Detect_ASAP',false);

%% Plot
ZS_threshold = 1;
Threshold = SynStats.amp_mean(3)+ZS_threshold*SynStats.amp_sigma(3);

[ax,h] = ROC_plot(detection_pk);
[ax,h] = ROC_plot(detection_pk,Threshold,'axes',ax,'h',h);
ROC_plot(detection_amp,Threshold);

% Normalized scale
xl = 100*detection_pk.tpr(1)/detection_pk.fpr(1);
xlim(ax,[0,xl]);
xlabel(ax,'False positive rate (Hz)/True peak rate (Hz)');
xt = 0:0.2:1;
xtl = cellfun(@num2str,num2cell(xt'),'UniformOutput',false);
xticks(ax,xt*xl);
xticklabels(ax,xtl);
% pause;

% Change back
xlim(ax,[0,100]);
xlabel(ax,'False positive rate (%)');
xticks(ax,'auto');
xticklabels(ax,'auto');

% Get conditional probability given detection
CProbGt = 100*interp1(detection_pk.ctrs_thr0,detection_pk.CProbEq',Threshold);
CProbEq = 100*interp1(detection_pk.ctrs_thr0,detection_pk.CProbEq',Threshold);
disp(num2str(CProbGt,['Given a detection greater than threshold: False: %.3g%%,', ...
    'Intermediate: %.3g%%, True: %.3g%%']));
disp(num2str(CProbEq,['Given a detection equal to threshold: False: %.3g%%,', ...
    'Intermediate: %.3g%%, True: %.3g%%']));
