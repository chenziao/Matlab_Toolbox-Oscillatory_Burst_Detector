close all;
clearvars -except input_directory

% load data
input_directory = '../example data';

if ~exist(input_directory,'dir'),   input_directory = pwd;    end
[file,input_directory] = uigetfile('','Select Input Data',input_directory);
if isequal(file,0), return; end
load(fullfile(input_directory,file));

save_directory = input_directory;
[~,file] = fileparts(file);
filename = [file,'_Charac'];

%% Input
% User input
% LFP_seg	% LFP segments data in a cell array
if exist('scale','var')
	LFP_seg = cellfun(@(x) scale*double(x),LFP_seg,'UniformOutput',false);
end
if ~exist('fs','var'),	fs = 1000;	end	% sampling freuency
disp(num2str(sum(cellfun(@length,LFP_seg))/fs/60,'Time length: %.1f minutes.'));

unit = 'μV';    % unit of LFP amplitude
freq_band = [30,100];  % signal frequency band (Hz)
fit_range = [20,250];	% frequency range for fitting a/f^b

%% Characterize PSD
% Parameter
db_threshold = 0.95;	% decibel threshold
autofit = true; % find range for fitting automatically (set to 2 to enable tightmode)
plotmode = 1;   % plot PSD fit results (set to 2 to plot autofit iterations)

[RealData,LFP_seg,flag] = getPSDestimation(LFP_seg,fs,'display',1);
if flag>-2
    [RealData,flag] = getPSDcomponents(RealData,freq_band,fit_range,...
        'tDB',db_threshold,'autofit',autofit,'plot',plotmode,'display',1);
end
switch flag
    case -2
        disp('The lengths of input LFP segments are too short for Fourier transform.');
    case -1
        disp('Frequency range for fitting is too narrow. Select wider range.');
    case 0
        disp('Power above the dB threshold not found. Suggestion: Lower the dB threshold.');
    case -3
        disp('Autofit failed: Signal occurs on the edge of auto fit range.');
end
if flag<=0,	return; end

%% Characterize Bursts
% Parameters
butter_order = 6;	% butterworth filter order
ZS_threshold = 2.0;	% amplitude Z-score threshold for burst detection
[RealData,stats] = getAllBurstAttr(RealData,LFP_seg,butter_order, ...
    ZS_threshold,'auto_threshold',0,'unit',unit,'fit_amp_gamma',1,'display',1);

%% Fit distributions
[RealData,stats] = getBurstDistributions(RealData,stats,'plot',1,'display',1,'plot_logscale',1);   % 'split_AP',1

%% Plot amplitude distribution
%%{
amp_ctrs = RealData.amp_ctrs;
gam_amp_sig = RealData.gam_amp_sig;
figure; hold on;
plot(amp_ctrs,raylpdf(amp_ctrs,RealData.pow_bg_filt^0.5),'r');
plot(amp_ctrs,gampdf(amp_ctrs,gam_amp_sig(1),gam_amp_sig(2)),'g');
plot(amp_ctrs,RealData.pdf_amp,'b');
plot(amp_ctrs,ncxpdf_cond_gam(amp_ctrs,gam_amp_sig(1),gam_amp_sig(2),RealData.pow_bg_filt^0.5),'m');
plot(amp_ctrs,gampdf(amp_ctrs,RealData.gam_amp_tot(1),RealData.gam_amp_tot(2)),'c');
legend({'Background','Fit signal','Composite','Fit composite','Gamma fit'});
xlabel(['Amplitude (',RealData.unit,')']);  ylabel('Probability density');
%}

%% Save results
[filename,save_directory] = uiputfile('*.mat','Save Results',fullfile(save_directory,filename));
if ~isequal(filename,0)
    save(fullfile(save_directory,filename),'-struct','RealData');
end
% Save burst statistics
%{
filename = [file,'_BurstStats'];
[filename,save_directory] = uiputfile('*.mat','Save Burst Statistics',fullfile(save_directory,filename));
if ~isequal(filename,0)
    save(fullfile(save_directory,filename),'-struct','stats');
end
%}
