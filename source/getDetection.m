function DetectionResult = getDetection(SynStats,BoundStats,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'SynStats',validStruct);
addRequired(p,'BoundStats',validStruct);
addOptional(p,'nthresh',50,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'Detect_AP',true,validTF);
parse(p,SynStats,BoundStats,varargin{:});

Custom_Detect = isfield(SynStats,'signal') && isfield(SynStats,'thresholds');
if Custom_Detect
    Detect_AP = false;
    target = SynStats.signal(:);    % signal being detected
    max_amp = max(target);  % maximum value of signal
    thresholds = SynStats.thresholds(:)';   % detection thresholds
    thresholds = sort(thresholds(thresholds<max_amp));
    nthresh = numel(thresholds);    % number of thresholds
    norm_factor = numel(target);    % normalizing factor for histogram count
    if isfield(SynStats,'index') && ~isempty(SynStats.index)
        idx = SynStats.index;   % time indices for detection
        F_idx = BoundStats.F_idx(idx);
        T_idx = BoundStats.T_idx(idx);
        if numel(F_idx)~=norm_factor
            error(['Input "signal" must have the same size as ' ...
                'the time points selected by "index".']);
        end
    else
        F_idx = BoundStats.F_idx;
        T_idx = BoundStats.T_idx;
        if numel(F_idx)~=norm_factor
            error('Input "signal" must have the same size as the synthetic data.');
        end
    end
    dur_prop = [sum(F_idx),0,sum(T_idx)]/norm_factor;
    dur_prop(2) = 1-sum(dur_prop);
else
    %% Select thresholds
    Detect_AP = p.Results.Detect_AP;    % whether detecting amplitude peaks
    nthresh = p.Results.nthresh;    % number of thresholds
    max_amp = SynStats.max_amp(3);  % maximum amplitude
    if Detect_AP
        gam_tot = SynStats.gam_AP_tot;
    else
        gam_tot = SynStats.gam_amp_tot;
    end
    thr_max = gaminv(1-1/nthresh^2,gam_tot(1),gam_tot(2));
    thr_max = min(thr_max,max_amp*(nthresh-2)/(nthresh-1));
    thresholds = linspace(0,thr_max,nthresh);   % detection thresholds
    %% Select target
    if Detect_AP
        target = SynStats.peaks{3};
        F_idx = BoundStats.F_idx(SynStats.tpks{3});
        T_idx = BoundStats.T_idx(SynStats.tpks{3});
        norm_factor = SynStats.syn_len;
    else
        target = SynStats.LFP_amp(:,3);
        F_idx = BoundStats.F_idx;
        T_idx = BoundStats.T_idx;
        norm_factor = SynStats.NT;
        dur_prop = BoundStats.dur_prop;
    end
end
%% Calculate FPR,TPR
% True/False/Intermediate classes
edge_thr = [thresholds,max_amp];
h_F = histcounts(target(F_idx),edge_thr);
h_T = histcounts(target(T_idx),edge_thr);
h_I = histcounts(target(~(F_idx|T_idx)),edge_thr);
h_fit = [h_F;h_I;h_T]/norm_factor;  % normalized histogram of F,I,T classes
fpr = fliplr(cumsum(fliplr(h_fit(1,:))));
tpr = fliplr(cumsum(fliplr(h_fit(3,:))));
if ~Detect_AP
    fpr = fpr/dur_prop(1);
    tpr = tpr/dur_prop(3);
end

%% Analysis
% For distribution stack plot
ctrs_thr = (edge_thr(1:end-1)+edge_thr(2:end))/2;
d_thr = diff(edge_thr);
pds = h_fit./d_thr; % probability density
cds = [zeros(size(ctrs_thr));cumsum(pds,1)]';	% cumulative density
% For conditional probability on detection
% conditioned on detection=threshold
CProbEq = [[1;0;0],pds./cds(:,4)',[0;0;1]];
ctrs_thr0 = [thresholds(1),ctrs_thr,max_amp];	% including two end points
rcdf = fliplr(cumtrapz(ctrs_thr0,fliplr([zeros(3,1),pds,zeros(3,1)]),2));	% reverse cdf
% conditioned on detection>threshold
CProbGt = [rcdf(:,1:end-1)./sum(rcdf(:,1:end-1),1),[0;0;1]];

% Youden's Index
[~,Youden] = max(tpr/tpr(1)-fpr/fpr(1));
if Detect_AP
    [~,Youden2] = max(tpr-fpr);
end
% TODO: AUC

%% Store result
varnames = {'Custom_Detect','Detect_AP','max_amp','nthresh','thresholds', ...
    'fpr','tpr','h_fit','ctrs_thr','cds','ctrs_thr0','CProbEq','CProbGt','Youden'};
if Detect_AP
    varnames = [varnames,'Youden2'];
end
DetectionResult = struct();
for i = 1:numel(varnames),	DetectionResult.(varnames{i})=eval(varnames{i});	end
DetectionResult.unit = SynStats.unit;

end
