close all;
clearvars -except working_directory;

if ~exist('working_directory','var')
    working_directory = fullfile('examples','example1');
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
% SynParam.fs = 1000;
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
[psd_bg,f,Nfft,fit_i] = SynParam.psd_background;
[pxx_syn,f] = pwelch(SynData.CompositeTrace,Nfft,[],Nfft,SynParam.fs);
pxx_syn([1,end]) = pxx_syn([1,end])*2;
[pxx_bg,f] = pwelch(SynData.BackgroundTrace,Nfft,[],Nfft,SynParam.fs);
pxx_bg([1,end]) = pxx_bg([1,end])*2;

figure(110);	hold on;
h(2) = plot(f,pow2db(pxx_bg),'b');
h(1) = plot(f,pow2db(psd_bg),'g');
legend(h,{'in vivo LFP of background','synthetic LFP of background'});
% axis([5,300,-10,40]); set(gca, 'XScale', 'log');
xlim([0,200]);
xlabel('Frequency (Hz)');   ylabel('PSD (dB/Hz)');

figure(111);	hold on;
plot(SynParam.rdat.f,pow2db(SynParam.rdat.PSD),'g');
plot(f,pow2db(pxx_syn),'b');
plot(f,pow2db(pxx_bg),'r');
xlim(SynParam.rdat.f(SynParam.rdat.i_fit_range));
% axis([30,120,2,20]);
legend({'in vivo LFP','synthetic LFP','synthetic backgound'});	%set(gca, 'XScale', 'log');
xlabel('Frequency (Hz)');   ylabel('PSD (dB/Hz)');

bump_f = SynParam.rdat.bump_f;
bump_i = find(f>bump_f(1) & f<bump_f(2));
intp = @(x,i) interp1(f,x,bump_f(i));
pow_psd = @(x) trapz([bump_f(1);f(bump_i);bump_f(2)],[intp(x,1);x(bump_i);intp(x,2)]);
pow_bg = pow_psd(pxx_bg);
pow_sig = pow_psd(pxx_syn)-pow_bg;
SNR = pow_sig/pow_bg;
%}

%% Characterize synthetic signal
%%{
SynParam.loadrealdata;
butter_order = SynParam.rdat.butter_order;
threshold = SynParam.rdat.threshold;
stop_perc = SynParam.rdat.stop_perc;
SynP = struct('fs',SynParam.fs,'nDB',SynParam.rdat.nDB,'bump_f',SynParam.rdat.bump_f, ...
    'T_length',SynParam.syn_len,'valid_seg_id',1,'f',f,'PSD_smoo',pxx_syn, ...
    'fit_i',fit_i,'fit_a',SynParam.rdat.fit_a,'fit_b',SynParam.rdat.fit_b);
[SynP,stats] = getAllBurstAttr(SynP,{SynData.CompositeTrace},butter_order, ...
    'threshold',threshold,'stop_perc',stop_perc,'display',1,'fit_ASamp',0);
SynP = getBurstDistributions(SynP,stats,'plot',1,'display',1,'plot_logscale',1);
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
xlabel('time (ms)');    ylabel('Amplitude (uV)');
legend(h,{'Synthetic signal','Synthetic LFP'},'Location','SouthEast');

cd(source_dir);
