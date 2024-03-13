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

%% Generate synthetic data
load(filename);
SynParam.syn_len = 5000;
% SynParam.bl_power_match = 1;
% SynParam.power_match = 1;
% SynParam.bg_pow_scale_whole = 1;
dist_type = SynParam.best_amp_dist_type;
if isempty(dist_type)
    return;
end
SynParam.amp_dist_type = dist_type;
% SynParam.amp_dist_type = 'exponential';
% SynParam.amp_dist_type = 'lognormal';
SynParam.amp_dist_type = 'gamma';
disp(['Amplitude distribution type being used: ',SynParam.amp_dist_type]);
SynParam.setup;
SynData = gen_synthetic(SynParam);

%% Save
[~,fname] = fileparts(filename);
fname = strrep(fname,'_SynParam','_SynData');
[filename,path] = uiputfile('*.mat','Save Results',fullfile(working_directory,fname));
if ~isequal(filename,0)
    working_directory = path;
    save(fullfile(working_directory,filename),'SynData','SynParam');
end

SynData.CompositeTrace = SynData.BackgroundTrace+SynData.SignalTrace;
SynParam.loadrealdata;

%% Plot PSD
%%{
[psd_bg,f,Nfft,bg_range_i] = SynParam.psd_background;
[pxx_syn,f] = pwelch(SynData.CompositeTrace,Nfft,[],Nfft,SynParam.fs);
pxx_syn([1,end]) = pxx_syn([1,end])*2;
[pxx_bg,f] = pwelch(SynData.BackgroundTrace,Nfft,[],Nfft,SynParam.fs);
pxx_bg([1,end]) = pxx_bg([1,end])*2;

figure(110);	hold on;
h(2) = plot(f,pow2db(pxx_bg),'b');
h(1) = plot(f,pow2db(psd_bg),'g');
legend(h,{'Observed background','Synthetic background'});
% axis([5,300,-10,40]); set(gca, 'XScale', 'log');
xlim([0,200]);
xlabel('Frequency (Hz)');   ylabel('PSD (dB/Hz)');

figure(111);	hold on;
plot(SynParam.rdat.f,pow2db(SynParam.rdat.PSD),'g');
plot(f,pow2db(pxx_syn),'b');
plot(f,pow2db(pxx_bg),'r');
xlim(SynParam.rdat.f(SynParam.rdat.i_fit_range));
set(gca, 'XScale', 'log');
xlim([10,400]);
legend({'Observed LFP','Synthetic LFP','Synthetic backgound'});
xlabel('Frequency (Hz)');   ylabel('PSD (dB/Hz)');

sig_range_f = SynParam.rdat.sig_range_f;
sig_range_i = find(f>sig_range_f(1) & f<sig_range_f(2));
intp = @(x,i) interp1(f,x,sig_range_f(i));
pow_psd = @(x) trapz([sig_range_f(1);f(sig_range_i);sig_range_f(2)], ...
    [intp(x,1);x(sig_range_i);intp(x,2)]);
pow_bg = pow_psd(pxx_bg);
pow_sig = pow_psd(pxx_syn)-pow_bg;
SNR = pow_sig/pow_bg;
disp(['Signal Noise Ratio: ',num2str(SNR,'%.2f '),' (',num2str(pow2db(SNR),'%.2f '),' dB)']);
%}

%% Characterize synthetic signal
%%{
SynParam.loadrealdata;
butter_order = SynParam.rdat.butter_order;
threshold = SynParam.rdat.threshold;
stop_perc = SynParam.rdat.stop_perc;
freq_resolution = SynParam.rdat.freq_resolution;
unit = SynParam.rdat.unit;
SynP = struct('fs',SynParam.fs,'tDB',SynParam.rdat.tDB,'sig_range_f',SynParam.rdat.sig_range_f, ...
    'T_length',SynParam.syn_len,'valid_seg_id',1,'f',f,'PSD_smoo',pxx_syn, ...
    'bg_range_i',bg_range_i,'fit_a',SynParam.rdat.fit_a,'fit_b',SynParam.rdat.fit_b);
if isfield(SynParam.rdat,'aperiodic_params')
    SynP.aperiodic_params = SynParam.rdat.aperiodic_params;
end
[SynP,stats] = getAllBurstAttr(SynP,{SynData.CompositeTrace},butter_order, ...
    'threshold',threshold,'stop_perc',stop_perc,'freq_resolution',freq_resolution, ...
    'unit',unit,'display',1,'fit_amp',0);
SynP = getBurstDistributions(SynP,stats,'plot',1,'plot_logscale',1,'display',1);
%}

%% Plot traces
tshow = [30,31];  % sec
ishow = round(tshow(1)*SynParam.fs):round(tshow(2)*SynParam.fs);
% tshow =1:5000;
figure; hold on;
% plot(ishow*SynParam.dt,SynData.BackgroundTrace(ishow),'Color',[1,.6,.8]);
h(2) = plot(ishow*SynParam.dt,SynData.CompositeTrace(ishow),'Color',0.5*[1,1,1]);
h(1) = plot(ishow*SynParam.dt,SynData.SignalTrace(ishow),'b');
axis tight;
xlabel('time (ms)');    ylabel(['Amplitude (',SynParam.rdat.unit,')']);
legend(h,{'Synthetic signal','Synthetic LFP'},'Location','SouthEast');

cd(source_dir);
