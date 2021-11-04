function DetectionResult = getDetection(SynStats,BoundStats,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'SynStats',validStruct);
addRequired(p,'BoundStats',validStruct);
addOptional(p,'nthresh',50,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParamValue(p,'Detect_ASAP',true,validTF);
parse(p,SynStats,BoundStats,varargin{:});

%% Select thresholds
Detect_ASAP = p.Results.Detect_ASAP;
nthresh = p.Results.nthresh;
max_amp = SynStats.max_amp(3);
if Detect_ASAP
    gamtot = SynStats.gam_pk_tot;
else
    gamtot = SynStats.gam_amp_tot;
end
thr_max = gaminv(1-1/nthresh^2,gamtot(1),gamtot(2));
thr_max = min(thr_max,max_amp*(nthresh-2)/(nthresh-1));
thresholds = linspace(0,thr_max,nthresh);

%% Calculate FPR,TPR
if Detect_ASAP
    target = SynStats.peaks{3};
    F_idx = BoundStats.F_idx(SynStats.tpks{3});
    T_idx = BoundStats.T_idx(SynStats.tpks{3});
else
    target = SynStats.LFP_amp(:,3);
    F_idx = BoundStats.F_idx;
    T_idx = BoundStats.T_idx;
end
% True/False/Intermediate classes
h_F = histc(target(F_idx),[thresholds,inf],1)';
h_T = histc(target(T_idx),[thresholds,inf],1)';
h_I = histc(target(~(F_idx|T_idx)),[thresholds,inf],1)';
h_fit = [h_F;h_I;h_T];	h_fit(:,end) = [];	% histogram of F,I,T classes
if Detect_ASAP
    h_fit = h_fit/SynStats.syn_len;
else
    h_fit = h_fit/SynStats.NT;
end
fpr = fliplr(cumsum(fliplr(h_fit(1,:))));
tpr = fliplr(cumsum(fliplr(h_fit(3,:))));
if ~Detect_ASAP
    fpr = fpr/BoundStats.dur_prop(1);
    tpr = tpr/BoundStats.dur_prop(3);
end

%% Analysis
% For distribution stack plot
edge_thr = [thresholds,max_amp];
ctrs_thr = (edge_thr(1:end-1)+edge_thr(2:end))/2;
d_thr = diff(edge_thr);
pds = bsxfun(@rdivide,h_fit,d_thr);	% probability density
cds = [zeros(size(ctrs_thr));cumsum(pds,1)]';	% cumulative density
% For conditional probability on detection
% conditioned on detection=threshold
CProbEq = [[1;0;0],bsxfun(@rdivide,pds,cds(:,4)'),[0;0;1]];
ctrs_thr0 = [0,ctrs_thr,max_amp];	% including two end points
rcdf = fliplr(cumtrapz(ctrs_thr0,fliplr([zeros(3,1),pds,zeros(3,1)]),2));	% reverse cdf
% conditioned on detection>threshold
CProbGt = [bsxfun(@rdivide,rcdf(:,1:end-1),sum(rcdf(:,1:end-1),1)),[0;0;1]];

% Youden's Index
[~,Youden] = max(tpr/tpr(1)-fpr/fpr(1));
if Detect_ASAP
    [~,Youden2] = max(tpr-fpr);
end
% TODO: AUC

%% Store result
varnames = {'Detect_ASAP','max_amp','nthresh','thresholds','fpr','tpr', ...
    'h_fit','ctrs_thr','cds','ctrs_thr0','CProbEq','CProbGt','Youden'};
if Detect_ASAP
    varnames = [varnames,'Youden2'];
end
DetectionResult = struct();
for i = 1:numel(varnames),	DetectionResult.(varnames{i})=eval(varnames{i});	end

end
