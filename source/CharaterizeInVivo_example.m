close all;
clearvars -except input_directory

% load data
input_directory = 'MultiElectrodeData';

if ~exist(input_directory,'dir'),   input_directory = pwd;    end
[file,input_directory] = uigetfile('','Select Input Data',input_directory);
if isequal(file,0), return; end
load(fullfile(input_directory,file));

save_directory = input_directory;
[~,file] = fileparts(file);
filename = [file,'_Charac'];

%% Input
% User input
% LFP_seg	% data in a cell array
if exist('scale','var')
	LFP_seg = cellfun(@(x) scale*double(x),LFP_seg,'UniformOutput',false);
end
if ~exist('fs','var'),	fs = 1000;	end	% sampling freuency
disp(num2str(sum(cellfun(@length,LFP_seg))/fs/60,'Time length: %.1f minutes.'));

freq_range = [30,100];  % frequency bound of interest (Hz)
fit_range = [20,250];	% range for fitting 1/f^a

%% Characterize PSD
% Parameter
outlier_threshold = 1;	% Decibel
autofit = true;

[RealData,LFP_seg,flag] = getPSDcomponents(LFP_seg,fs,freq_range,fit_range,...
    'nDB',outlier_threshold,'autofit',autofit,'plot',1,'display',1);
switch flag
    case -2
        disp('Input LFP data is too short. Segment durations are expected to be longer.');
    case -1
        disp('Frequency range for fitting is too narrow. Select wider range.');
    case 0
        disp('Suggestion: Lower the outlier threshold.');
    case -3
        disp('Autofit failed: Signal occurs on the edge of auto fit range.');
end
if flag<=0,	return; end

%% Characterize Bursts
% Parameters
ZS_threshold = 2.0;	% AS amp Z-score threshold for burst detection
butter_order = 6;	% filter order
[RealData,stats] = getAllBurstAttr(RealData,LFP_seg,butter_order,ZS_threshold,'display',1);

%% Fit distributions
[RealData,stats] = getBurstDistributions(RealData,stats,'plot',1,'display',1,'plot_logscale',1);    % 'with_lo',1

%% Plot amplitude distribution
%%{
as_amp = RealData.as_amp;
gam_amp_sig = RealData.gam_amp_sig;
figure; hold on;
plot(as_amp,raylpdf(as_amp,RealData.pow_bg_filt^0.5),'r');
plot(as_amp,gampdf(as_amp,gam_amp_sig(1),gam_amp_sig(2)),'g');
plot(as_amp,RealData.pdf_amp,'b');
plot(as_amp,ncxpdf_cond_gam(as_amp,gam_amp_sig(1),gam_amp_sig(2),RealData.pow_bg_filt^0.5),'m');
plot(as_amp,gampdf(as_amp,RealData.gam_amp_tot(1),RealData.gam_amp_tot(2)),'c');
legend({'Background','Fit Signal','Composite','Fit Composite','Fit Gamma'});
xlabel('AS amplitude (\muV)');	ylabel('Probability density');
%}

%% Save results
[filename,save_directory] = uiputfile('*.mat','Save Results',fullfile(save_directory,filename));
if ~isequal(filename,0)
    save(fullfile(save_directory,filename),'-struct','RealData');
end
