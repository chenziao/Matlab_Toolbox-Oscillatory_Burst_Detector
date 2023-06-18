function [axes,h] = ROC_plot(DetectionResult,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'DetectionResult',validStruct);
addOptional(p,'Threshold',[],@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
addParameter(p,'plot_dist',true,validTF);
addParameter(p,'axes',[]);
addParameter(p,'h',[]);
parse(p,DetectionResult,varargin{:});

Threshold = p.Results.Threshold;
plot_dist = logical(p.Results.plot_dist);
axes = p.Results.axes;
h = p.Results.h;
newplot = isempty(Threshold) || isempty(axes) || isempty(h);
if isempty(axes)
    figure;	ax = gca;
    axes = ax;
    if plot_dist
        figure;	axh = gca;
        axes = [ax,axh];
    end
else
    ax = axes(1);
    if plot_dist,   axh = axes(2);  end
end
if isempty(h)
    h = cell(plot_dist+1,2);
    [h{:,2}] = deal(0);
elseif size(h,1)<plot_dist+1
    h{2,2} = 0;
end
if isempty(Threshold)
    Threshold = DetectionResult.max_amp;
else
    Threshold = min(Threshold,DetectionResult.max_amp);
    for i = 1:plot_dist+1
        delete(h{i,2}(h{i,2}~=0 & ishandle(h{i,2})));
    end
end

Custom_Detect = isfield(DetectionResult,'Custom_Detect') && DetectionResult.Custom_Detect;
Detect_AP = DetectionResult.Detect_AP;
fpr = DetectionResult.fpr;
tpr = DetectionResult.tpr;

%% Plot ROC
hold(ax,'on');
if newplot
    pfpr = 100/fpr(1)*fpr;
    ptpr = 100/tpr(1)*tpr;
    if Detect_AP, colororder(ax,{'k','k'});	end
    hh = zeros(1,3);
    hh(1) = plot(ax,pfpr,ptpr,'b');
    hh(2) = plot(ax,[0,100],[0,100],':','color',0.5*[1,1,1]);
    YD = DetectionResult.Youden;
    plot(ax,pfpr(YD)*[1,1],[pfpr(YD),ptpr(YD)],'b:');
    hh(3) = plot(ax,pfpr(YD),ptpr(YD),'bd');
    xlim(ax,[0,100]);   ylim(ax,[0,100]);
    xlabel(ax,'False positive rate (%)');
    ylabel(ax,'True positive rate (%)');
    if Detect_AP
        hh = [hh(1:2),0,hh(3),0];
        yyaxis(ax,'right')
        hh(3) = plot(ax,[0,100*tpr(1)/fpr(1)],[0,tpr(1)],'--','color',0.5*[1,1,1]);
        YD2 = DetectionResult.Youden2;
        plot(ax,pfpr(YD2)*[1,1],[fpr(YD2),tpr(YD2)],'b--');
        hh(5) = plot(ax,pfpr(YD2),tpr(YD2),'bs');
        ylim(ax,[0,tpr(1)]);
        ylabel(ax,'True positive rate (Hz)');
        yyaxis(ax,'left');
    end
    h{1,1} = hh;
end

FPR = 100/fpr(1)*interp1(DetectionResult.thresholds,fpr,Threshold,[],0);
TPR = 100/tpr(1)*interp1(DetectionResult.thresholds,tpr,Threshold,[],0);
h{1,2} = plot(ax,FPR,TPR,'ro');
hh = [h{1,1},h{1,2}];

if Detect_AP
    leg = {'ROC','Random','TP rate = FP rate (Hz)','Maximum Youden''s index', ...
        'Youden''s index in Hz','Operating point'};
else
    leg = {'ROC','Random','Maximum Youden''s index','Operating point'};
end
legend(ax,hh,leg,'Location','SouthEast');

%% Distribution stack plot
if plot_dist
    % Define patches
    nthresh = DetectionResult.nthresh;
    ctrs_thr = DetectionResult.ctrs_thr;
    v2 = [[repmat(ctrs_thr',4,1),DetectionResult.cds(:)]; ...
        [DetectionResult.thresholds(1),0;DetectionResult.max_amp,0]];
    f2 = zeros(3,nthresh*2+2);
    endpt = 4*nthresh+[1,2];
    for i = 1:3
        f2(i,:) = [endpt(1),i*nthresh+1:(i+1)*nthresh,endpt(2),i*nthresh:-1:(i-1)*nthresh+1];
    end
    CDS = @(i,j) v2(f2(i,1:nthresh+2),j);
    % Plot distribution
    cmap = [1,0,0;.6,.6,.6;0,0,1];
    cmap1 = @(i,a) a*cmap(i,:)+(1-a)*[1,1,1];
    hold(axh,'on');
    if newplot
        % plot face
        hh = zeros(1,3);
        for i = 1:3
            hh(i) = patch(axh,'Faces',f2(i,:),'Vertices',v2,'FaceColor',cmap1(i,0.4),'EdgeColor','none');
        end
        h{2,1} = hh;
        % plot lines
        cmap2 = cmap([1,3,3],:);
        for i = 3:-1:1
            plot(axh,CDS(i,1),CDS(i,2),'color',cmap2(i,:));
        end
        unit = DetectionResult.unit;
        if Detect_AP
            xlabel(axh,['Amplitude peak (',unit,')']);
            ylabel(axh,['Rate density (Hz/',unit,')']);
        else
            if Custom_Detect
                xlabel(axh,'Custom detection signal');
            else
                xlabel(axh,['Amplitude (',unit,')']);
            end
            ylabel(axh,['Probability density (',unit,'^{-1})']);
        end
    end
    
    % plot threshold
    dsthr = arrayfun(@(i) interp1(CDS(i,1),CDS(i,2),Threshold,[],0),1:3);
    v2 = [v2;[repmat(Threshold,6,1),[0,dsthr(1:2),dsthr]']];
    i_thr = find([ctrs_thr,Inf]>Threshold,1);
    thrpt = reshape(4*nthresh+2+(1:6),3,2);
    f22 = zeros(3,(nthresh-i_thr+1)*2+3);
    for i = 1:3
        f22(i,:) = [thrpt(i,:),i*nthresh+i_thr:(i+1)*nthresh,endpt(2),i*nthresh:-1:(i-1)*nthresh+i_thr];
    end
    % plot face
    h{2,2} = zeros(1,4);
    for i = 1:3
        h{2,2}(i+1) = patch(axh,'Faces',f22(i,:),'Vertices',v2,'FaceColor',cmap1(i,0.7),'EdgeColor','none');
    end
    h{2,2}(1) = plot(axh,Threshold*[1,1],[0,dsthr(3)],'y');
    hh = [h{2,1},h{2,2}(1)];
    % bring lines to top
    hhh = get(axh,'Children');
    idx_lines = strcmp(get(hhh,'Type'),'line') & hhh~=hh(end);
    set(axh,'Children',[hhh(idx_lines);hhh(~idx_lines)]);
    
    leg = {'False peaks','Intermediate peaks','True peaks','Detection threshold'};
    legend(axh,hh,leg,'Location','NorthEast');
end

end
