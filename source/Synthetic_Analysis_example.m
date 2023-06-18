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
disp(num2str(SynParam.syn_len,'Synthetic LFP length: %.1f seconds'));

%% Primary Process
[SynStats,SynData] = getSyntheticStats(SynParam,SynData); % Update SynData with field "CompositeTrace"

%% Preselected bounds
BoundStats = getPreselectedBounds(SynParam,SynStats);
disp(num2str(SynParam.rdat.tDB,'PSD decibel threshold: %.2f dB'));
disp(num2str(BoundStats.auto_thetaU,['Auto upper bound: %.4g ',SynStats.unit]));
disp(num2str(BoundStats.pre_bound,['Preselected lower/upper bound: %.4g, %.4g ',SynStats.unit]));
disp(num2str(BoundStats.ZS_bound,'Z-score of background amplitude: %.2f, %.2f'));
disp(num2str(BoundStats.nsigma_bound,'¦È_L = %.2f¦Ò, ¦È_U = %.2f¦Ò'));

%% Detection
detection_AP = getDetection(SynStats,BoundStats,'Detect_AP',true);
detection_amp = getDetection(SynStats,BoundStats,'Detect_AP',false);
% % Example for customized detection method
% LFP_filt = filtfilt(SynParam.sos,SynParam.g,SynData.CompositeTrace);    % filtered trace
% [trs,ttr] = findpeaks(-LFP_filt);  % detect on amplitude only at troughs
% Custom_Detect = struct('signal',trs,'thresholds',detection_amp.thresholds); % specify signal and thresholds
% Custom_Detect.index = ttr; % time indices are required for identifying ground truth labels
% detection_amp = getDetection(Custom_Detect,BoundStats);

%% Plot
ZS_threshold = 1;
Threshold = SynStats.amp_mean(3)+ZS_threshold*SynStats.amp_sigma(3);

[ax,h] = ROC_plot(detection_AP);
[ax,h] = ROC_plot(detection_AP,Threshold,'axes',ax,'h',h);
ROC_plot(detection_amp,Threshold);

% Normalized FPR scale
%{
xl = 100*detection_AP.tpr(1)/detection_AP.fpr(1);
xlim(ax(1),[0,xl]);
xlabel(ax(1),'False positive rate (Hz)/True peak rate (Hz)');
xt = 0:0.2:1;
xtl = cellfun(@num2str,num2cell(xt'),'UniformOutput',false);
xticks(ax(1),xt*xl);
xticklabels(ax(1),xtl);
pause;
% Change back
xlim(ax(1),[0,100]);
xlabel(ax(1),'False positive rate (%)');
xticks(ax(1),'auto');
xticklabels(ax(1),'auto');
%}

% Get conditional probability given detection
CProbGt = 100*interp1(detection_AP.ctrs_thr0,detection_AP.CProbEq',Threshold);
CProbEq = 100*interp1(detection_AP.ctrs_thr0,detection_AP.CProbEq',Threshold);
disp(num2str(CProbGt,['Given a detection greater than threshold: False: %.3g%%, ', ...
    'Intermediate: %.3g%%, True: %.3g%%']));
disp(num2str(CProbEq,['Given a detection equal to threshold: False: %.3g%%, ', ...
    'Intermediate: %.3g%%, True: %.3g%%']));

%% Save results
SynResult = struct();
varnames = {'unit','syn_len','NT','npk','AP_rates','sigma_syn','max_amp','amp_mean', ...
    'amp_sigma','amp_skew','amp_kurt','gam_amp_syn','gam_AP_syn','gam_amp_tot','gam_AP_tot'};
for i = 1:numel(varnames)
    SynResult.(varnames{i}) = SynStats.(varnames{i});
end
varnames = {'pre_bound','nsigma_bound','ZS_bound','auto_thetaU','dur_prop','APrate_fit','APrate_sig_fit'};
for i = 1:numel(varnames)
    SynResult.(varnames{i}) = BoundStats.(varnames{i});
end
SynResult.detection_AP = detection_AP;
SynResult.detection_amp = detection_amp;

[~,fname] = fileparts(filename);
fname = strrep(strrep(fname,'_SynParam','_SynData'),'_SynData','_SynResult');
[filename,path] = uiputfile('*.mat','Save Results',fullfile(working_directory,fname));
if ~isequal(filename,0)
    working_directory = path;
    save(fullfile(working_directory,filename),'SynResult','SynParam');
end
