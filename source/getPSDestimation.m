function [results,LFP_seg,flag] = getPSDestimation(LFP_seg,fs,varargin)
p = inputParser;
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'LFP_seg',@(x) validateattributes(x,{'cell','numeric'},{'vector'}));
addRequired(p,'fs',validPositiveScalar);
% window size estimating PSD (NFFT)
addParameter(p,'nfft',8192,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
% window size for smoothing PSD (Hz)
addParameter(p,'smooth_win',2,validPositiveScalar);
addParameter(p,'display',false,validTF);
addParameter(p,'msgfcn',[]);
parse(p,LFP_seg,fs,varargin{:});

nfft = max(p.Results.nfft,32);  % at least 32
noverlap = floor(nfft/2);   % 50% overlap
[~,f] = pwelch(zeros(nfft,1),nfft,noverlap,nfft,fs);

%% Preprocess parameters
flag = -2;
results = [];
if ~iscell(LFP_seg)
    LFP_seg = {LFP_seg};
end
for i = 1:numel(LFP_seg)
    if ~isvector(LFP_seg{i}) || ~isnumeric(LFP_seg{i})
        error('Expect input data to be a cell vector, each element in which is a numeric vector.');
    end
    if ~isreal(LFP_seg{i})
        error('Expect input data to be real valued.');
    end
    LFP_seg{i} = double(LFP_seg{i}(:));
end
LFP_seg = LFP_seg(:);
LFP_mean = sum(cellfun(@sum,LFP_seg))/sum(cellfun(@length,LFP_seg));
LFP_seg = cellfun(@(x) x-LFP_mean,LFP_seg,'UniformOutput',0);

%% remove short segments
size_seg = reshape(cellfun(@length,LFP_seg),[],1);
validseg = size_seg>=nfft;
T_length = sum(size_seg(validseg))/fs;
valid_seg_id = int32(find(validseg));
if isempty(valid_seg_id),	return;	end
flag = -1;

%% calculate PSD
nwindows = floor((size_seg(validseg)-noverlap)/(nfft-noverlap));	% # of windows
pxxs = cell2mat(reshape(cellfun(@(x) pwelch(x,nfft,noverlap,nfft,fs),LFP_seg(validseg),'UniformOutput',0),1,[]));
PSD = pxxs*nwindows/sum(nwindows);
PSD(1) = 0;

% smooth PSD
span = find(f<p.Results.smooth_win,1,'last');
PSD_smoo = movmean(PSD,span);

%% Results
varnames = {'T_length','fs','valid_seg_id','f','PSD','PSD_smoo'};
results = eval(['struct(',strjoin(cellfun(@(x) ['''',x,''',',x],varnames,'UniformOutput',0),','),')']);

%% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    dispfcn(num2str([length(valid_seg_id),length(validseg)], ...
        '%d out of %d segments are long enough to use.'));
end

end