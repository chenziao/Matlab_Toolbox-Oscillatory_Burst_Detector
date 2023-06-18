function [results,stats] = getBurstDistributions(results,stats,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',validStruct);
addRequired(p,'stats',validStruct);
% number of bins for histograms
addOptional(p,'nbins',30,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'plot',false,validTF);
% plot variables distribution in logscale
addParameter(p,'plot_logscale',false,validTF);
% plot multiple fits of amplitude peak distribution
addParameter(p,'split_AP',false,validTF);
addParameter(p,'display',false,validTF);
addParameter(p,'axes',[]);
addParameter(p,'msgfcn',[]);
parse(p,results,stats,varargin{:});

peaks = stats.peaks;
npeak = stats.npeak;
nvalid = stats.nvalid;
threshold = results.threshold;
unit = results.unit;
nbins = p.Results.nbins;
split_AP = p.Results.split_AP;

%% AP distribution
mean_AP = mean(peaks);
% gamma
gam_AP = gamfit(peaks);
% exponential
exp_AP = expfit(stats.AP-threshold);
% lognormal
log_AP = log(peaks);
[logn_AP(1),logn_AP(2)] = normfit(log_AP);
% histogram
threshold_dist = max(results.mean_amp,min(peaks));	% threshold to split AP into two parts for fitting
[hist_AP,AP_edges] = histcounts(peaks,2*nbins);
prop_pks = cumsum(hist_AP); prop_pks = prop_pks/prop_pks(end);
ind1 = find(prop_pks<=1-1e-5,1,'last');
ind2 = find(medfilt1(hist_AP,5),1,'last');
switch sum([isempty(ind1),isempty(ind2)])
    case 0
        maxpeak = AP_edges(min(ind1,ind2)+1);
    case 1
        maxpeak = AP_edges([ind1,ind2]+1);
    case 2
        maxpeak = max(peaks);
end
maxpeak = min(maxpeak,10*results.mean_amp);
AP_edges = linspace(0,maxpeak,2*nbins+1);
AP_ctrs = bincenters(AP_edges);
hist_AP = histcounts(peaks,AP_edges);
pdf_AP = hist_AP/npeak/binwidth(AP_edges);
rd_AP = hist_AP/results.T_length/binwidth(AP_edges);

%% CN distribution
mean_CN = mean(stats.CN);
% lognormal
log_CN = log(stats.CN);
[logn_CN(1),logn_CN(2)] = normfit(log_CN);
% histogram
[pdf_CN,CN_edges] = histcounts(stats.CN,nbins,'Normalization','pdf');
log_CN_edges = linspace(min(log_CN),max(log_CN),nbins);
pdf_log_CN = histcounts(log_CN,log_CN_edges,'Normalization','pdf');
% empirical cdf of CN
pdf_log_CN_nonzero = max(pdf_log_CN,1/nvalid/binwidth(log_CN_edges));
cdf_log_CN = cumsum(pdf_log_CN_nonzero);
cdf_log_CN = [0,cdf_log_CN/cdf_log_CN(end)];
log_CNs = log_CN_edges;

%% BF distribution
mean_BF = mean(stats.BF);
% lognormal
log_BF = log(stats.BF);
logn_BF = fitdist_trunc(log_BF,@normpdf,@normcdf,log(results.sig_range_f),[mean(log_BF),std(log_BF)]);
% histogram
rng(0);
burst_freq_ctn = stats.BF+(rand(size(stats.BF))-0.5)*stats.df;  % dequantize to continuous value
[pdf_BF,BF_edges] = histcounts(burst_freq_ctn,nbins,'Normalization','pdf');
log_BF_edges = linspace(log(min(stats.BF)-stats.df/2),log(max(stats.BF)+stats.df/2),nbins);
pdf_log_BF = histcounts(log(burst_freq_ctn),log_BF_edges,'Normalization','pdf');
% empirical cdf of BF
BF_bw = binwidth(log_BF_edges);
pdf_scale = diff(normcdf(log_BF_edges([1,end]),logn_BF(1),logn_BF(2)));
pdf_scale = pdf_scale*nvalid/(nvalid+sum(pdf_log_BF==0));
pdf_log_BF_nonzero = pdf_scale*max(pdf_log_BF,1/nvalid/BF_bw);
freq_edge = logn_BF(1)+4*logn_BF(2)*[-1,1];
log_BF_l = fliplr(log_BF_edges(1)-BF_bw:-BF_bw:freq_edge(1));
log_BF_r = log_BF_edges(end)+BF_bw:BF_bw:freq_edge(2);
pdf_log_BF_l = normpdf(log_BF_l+BF_bw/2,logn_BF(1),logn_BF(2));
pdf_log_BF_r = normpdf(log_BF_r-BF_bw/2,logn_BF(1),logn_BF(2));
pdf_log_BF_wide = [pdf_log_BF_l,pdf_log_BF_nonzero,pdf_log_BF_r];
cdf_log_BF = cumsum(pdf_log_BF_wide);
cdf_log_BF = [0,cdf_log_BF/cdf_log_BF(end)];
log_BFs = [log_BF_l,log_BF_edges,log_BF_r];

%% caculate covariance
X = [stats.AP,stats.CN,stats.BF];
CorrLin = corrcoef(X);
CorrLog = corrcoef(log(X));

%% Add results
varnames = {'mean_AP','gam_AP','exp_AP','logn_AP','mean_CN','logn_CN','mean_BF','logn_BF', ...
    'CN_edges','pdf_CN','log_CNs','cdf_log_CN','BF_edges','pdf_BF','log_BFs','cdf_log_BF', ...
    'AP_edges','AP_ctrs','hist_AP','rd_AP','pdf_AP','CorrLin','CorrLog'};
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

%% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    dispfcn(num2str([mean_AP,threshold,mean(stats.AP)], ...
        ['Mean amplitude peak: %.5g ',unit,', amplitude threshold: %.5g ', ...
        unit,', mean valid burst amplitude peak: %.5g ',unit]));
    dispfcn([num2str(exp_AP,['Exponential fit mean: %.5g ',unit]), ...
        num2str(logn_AP,', Lognormal fit mean/(std): %.3g/(%.3g)'), ...
        num2str(gam_AP,[', Gamma fit (k,',char(952),'): (%.3f,%.5g)'])]);
    dispfcn(num2str([mean_CN,logn_CN],...
        'Mean cycle number: %.2f, Lognormal fit mean/(std): %.3g/(%.3g)'));
    dispfcn(num2str([mean_BF,logn_BF],...
        'Mean burst frequency: %.2f Hz, Lognormal fit mean/(std): %.3g/(%.3g)'));
    dispfcn('Correlation coefficients: ');
    dispfcn(num2str(CorrLin,'%.3f '));
    dispfcn('Logscale correlation coefficients: ');
    dispfcn(num2str(CorrLog,'%.3f '));
end
if p.Results.plot
    ax = p.Results.axes;    axflag = isempty(ax);
    plot_logscale = p.Results.plot_logscale;
    % AP
    % separate two parts
    if split_AP
        nbins_lo = floor(2*nbins*threshold_dist/maxpeak);
        split_AP = nbins_lo>=2;
    end
    % linear scale
    if split_AP
        i_AP_hi = peaks>=threshold_dist;
        logn_AP_hi = fitdist_trunc(log_AP,@normpdf,@normcdf,[log(threshold_dist),Inf], ...
            [mean(log_AP(i_AP_hi)),std(log_AP(i_AP_hi))]);
        npk = sum(i_AP_hi);
        AP_bw = threshold_dist/nbins_lo;
        nbins_hi = ceil((maxpeak-threshold_dist)/AP_bw);
        AP_edges = linspace(0,threshold_dist+nbins_hi*AP_bw,nbins_lo+nbins_hi+1);
        AP_lo = AP_edges(1:nbins_lo+1);
        AP_hi = AP_edges(nbins_lo+1:end);
        pdf_AP = histcounts(peaks,AP_edges,'Normalization','pdf');
        pdf_AP_lo = pdf_AP(1:nbins_lo);
        pdf_AP_hi = pdf_AP(nbins_lo+1:end);
    end
    threshold_idx = find(AP_edges>=threshold,1);
    % log scale
    if plot_logscale
        log_threshold = log(threshold_dist);
        log_AP_bw = (log(maxpeak)-log_threshold)/nbins;
        nbins_lo = ceil((log_threshold-log(min(peaks)))/log_AP_bw);
        log_AP_edges = linspace(log_threshold-nbins_lo*log_AP_bw,log(maxpeak),nbins_lo+nbins+1);
        pdf_log_AP = histcounts(log_AP,log_AP_edges,'Normalization','pdf');
        log_AP_ctrs = bincenters(log_AP_edges);
        log_AP_exp = exp(log_AP_ctrs);
        if split_AP
            log_AP_lo = log_AP_edges(1:nbins_lo+1);
            log_AP_hi = log_AP_edges(nbins_lo+1:end);
            pdf_log_AP_lo = pdf_log_AP(1:nbins_lo);
            pdf_log_AP_hi = pdf_log_AP(nbins_lo+1:end);
        end
        threshold_log_idx = find(log_AP_edges>=log(threshold),1);
    end
    if axflag
        figure;
        if plot_logscale,	subplot(211);   end
        ax(1) = gca;
    end
    hold(ax(1),'on');
    if split_AP
        histogram(ax(1),'BinEdges',AP_hi,'BinCounts',pdf_AP_hi);
        histogram(ax(1),'BinEdges',AP_lo,'BinCounts',pdf_AP_lo,'FaceColor','b');
        sig_lo = fitdist_trunc(peaks,@maxwboltzpdf,@maxwboltzcdf, ...
            [0,threshold_dist],(pi/8)^0.5*mean(peaks(~i_AP_hi)));
        plot(ax(1),AP_ctrs,(1-npk/npeak)/maxwboltzcdf(threshold_dist,sig_lo)*maxwboltzpdf(AP_ctrs,sig_lo),'c');
        plot(ax(1),AP_ctrs(threshold_idx:end),stats.ndetect/npeak ...
            *exppdf(AP_ctrs(threshold_idx:end),exp_AP)/(1-expcdf(threshold,exp_AP)),'m','LineWidth',2);
        plot(ax(1),AP_ctrs,npk/npeak/(1-logncdf(threshold_dist,logn_AP_hi(1),logn_AP_hi(2))) ...
            *lognpdf(AP_ctrs,logn_AP_hi(1),logn_AP_hi(2)),'y');
    else
        histogram(ax(1),'BinEdges',AP_edges,'BinCounts',pdf_AP);
    end
    plot(ax(1),AP_ctrs,gampdf(AP_ctrs,gam_AP(1),gam_AP(2)),'r');
    xlim(ax(1),AP_edges([1,end]));
    xlabel(ax(1),['Amplitude peak (',unit,')']);
    ylabel(ax(1),'Probability density');
    if split_AP
        leg = legend(ax(1),{'High amplitude','Low amplitude',...
            'Maxwell-Boltzmann fit','Exponential fit',...
            'Lognormal fit','Gamma fit'},'Location','NorthEast');
    else
        leg = legend(ax(1),{'Observed','Gamma fit'},'Location','NorthEast');
    end
    set(leg,'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(4) = gca;	end
        hold(ax(4),'on');
        if split_AP
            histogram(ax(4),'BinEdges',log_AP_hi,'BinCounts',pdf_log_AP_hi);
            histogram(ax(4),'BinEdges',log_AP_lo,'BinCounts',pdf_log_AP_lo,'FaceColor','b');
            plot(ax(4),log_AP_ctrs,(1-npk/npeak)/maxwboltzcdf(threshold_dist,sig_lo) ...
                *log_AP_exp.*maxwboltzpdf(log_AP_exp,sig_lo),'c');
            plot(ax(4),log_AP_ctrs,npk/npeak/(1-logncdf(threshold_dist,logn_AP_hi(1),logn_AP_hi(2))) ...
                *normpdf(log_AP_ctrs,logn_AP_hi(1),logn_AP_hi(2)),'y');
            plot(ax(4),log_AP_ctrs(threshold_log_idx:end),stats.ndetect/npeak ...
                *log_AP_exp(threshold_log_idx:end).*exppdf(log_AP_exp(threshold_log_idx:end),exp_AP) ...
                /(1-expcdf(threshold,exp_AP)),'m','LineWidth',2);
        else
            histogram(ax(4),'BinEdges',log_AP_edges,'BinCounts',pdf_log_AP);
        end
        plot(ax(4),log_AP_ctrs,log_AP_exp.*gampdf(log_AP_exp,gam_AP(1),gam_AP(2)),'r');
        xlim(ax(4),log_AP_edges([1,end]));
        xlabel(ax(4),'Amplitude peak (log scale)');
        ylabel(ax(4),'Probability density');
    end
    % CN
    if axflag
        figure;
        if plot_logscale,	subplot(211);   end
        ax(2) = gca;
    end
    hold(ax(2),'on');
    histogram(ax(2),'BinEdges',CN_edges,'BinCounts',pdf_CN);
    CN_ctrs = bincenters(CN_edges);
    plot(ax(2),CN_ctrs,lognpdf(CN_ctrs,logn_CN(1),logn_CN(2)),'r');
    xlim(ax(2),CN_edges([1,end]));
    xlabel(ax(2),'Cycle number');
    ylabel(ax(2),'Probability density');
    set(legend(ax(2),{'Observed','Lognormal fit'}),'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(5) = gca;	end
        hold(ax(5),'on');
        histogram(ax(5),'BinEdges',log_CN_edges,'BinCounts',pdf_log_CN);
        log_CN_ctrs = bincenters(log_CN_edges);
        plot(ax(5),log_CN_ctrs,normpdf(log_CN_ctrs,logn_CN(1),logn_CN(2)),'r');
        xlim(ax(5),log_CNs([1,end]));
        xlabel(ax(5),'Cycle number (log scale)');
        ylabel(ax(5),'Probability density');
    end
    % BF
    if axflag
        figure;
        if plot_logscale,	subplot(211);   end
        ax(3) = gca;
    end
    hold(ax(3),'on');
    histogram(ax(3),'BinEdges',BF_edges,'BinCounts',pdf_BF);
    BF_ctrs = bincenters(BF_edges);
    plot(ax(3),BF_ctrs,lognpdf(BF_ctrs,logn_BF(1),logn_BF(2)),'r');
    xlim(ax(3),BF_edges([1,end]));
    xlabel(ax(3),'Burst frequency (Hz)');
    ylabel(ax(3),'Probability density');
    set(legend(ax(3),{'Observed','Lognormal fit'}),'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(6) = gca;	end
        hold(ax(6),'on');
        histogram(ax(6),'BinEdges',log_BF_edges,'BinCounts',pdf_log_BF);
        log_BF_ctrs = bincenters(log_BF_edges);
        plot(ax(6),log_BF_ctrs,normpdf(log_BF_ctrs,logn_BF(1),logn_BF(2)),'r');
        h = plot(ax(6),bincenters(log_BFs),pdf_log_BF_wide,'c');
        xlim(ax(6),log_BFs([1,end]));
        xlabel(ax(6),'Burst frequency (log scale)');
        ylabel(ax(6),'Probability density');
        set(legend(ax(6),h,'Empirical + fit'),'color','none');
    end
end

end
