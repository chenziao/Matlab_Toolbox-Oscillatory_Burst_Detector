function [results,stats] = getBurstDistributions(results,stats,varargin)
p = inputParser;
validStruct = @(x) validateattributes(x,{'struct'},{'scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'results',validStruct);
addRequired(p,'stats',validStruct);
% number of bins for histograms
addOptional(p,'nbins',30,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParamValue(p,'plot',false,validTF);
% plot variables distribution in logscale
addParamValue(p,'plot_logscale',false,validTF);
% plot multiple types of amplitude distribution
addParamValue(p,'with_lo',false,validTF);
addParamValue(p,'display',false,validTF);
addParamValue(p,'axes',[]);
addParamValue(p,'msgfcn',[]);
parse(p,results,stats,varargin{:});

peaks = stats.peaks;
npeak = stats.npeak;
nvalid = stats.nvalid;
threshold = results.threshold;
nbins = p.Results.nbins;
with_lo = p.Results.with_lo;

% AP distribution
meanAP = mean(peaks);
% gamma
gamAP = gamfit(peaks);
% exponential
expAP = expfit(stats.AP-threshold);
% lognormal
log_AP = log(peaks);
[logAP(1),logAP(2)] = normfit(log_AP);
% histogram
threshold_dist = max(results.mean_amp,min(peaks));	% AP threshold to separate two parts for fitting
[hist_pks,pk_amp] = hist(peaks,2*nbins);
prop_pks = cumsum(hist_pks);	prop_pks = prop_pks/prop_pks(end);
ind1 = find(prop_pks<=1-1e-5,1,'last');
ind2 = find(medfilt1(hist_pks,5),1,'last');
switch sum([isempty(ind1),isempty(ind2)])
    case 0
        maxpeak = pk_amp(min(ind1,ind2))+binwidth(pk_amp)/2;
    case 1
        maxpeak = pk_amp([ind1,ind2])+binwidth(pk_amp)/2;
    case 2
        maxpeak = max(peaks);
end
maxpeak = min(maxpeak,10*results.mean_amp);
pk_amp = linspace(0,maxpeak,2*nbins+1);
pk_amp = (pk_amp(1:end-1)+pk_amp(2:end))/2;
hist_pks = hist(peaks,pk_amp);
pdf_pks = hist_pks/npeak/binwidth(pk_amp);
rd_pks = hist_pks/results.T_length/binwidth(pk_amp);

% CN distribution
meanCN = mean(stats.CN);
% lognormal
log_CN = log(stats.CN);
[logCN(1),logCN(2)] = normfit(log_CN);
% histogram
[hist_cyc,cyc] = hist(stats.CN,nbins);
log_cyc = linspace(min(log_CN),max(log_CN),nbins);
hist_cyc_log = hist(log_CN,log_cyc);
idx = hist_cyc_log==0;	hist_cyc_log(idx) = 1;
hist_cyc_log = nvalid/(nvalid+sum(idx))*hist_cyc_log;
% empirical cdf of CN
cdf_cyc = cumtrapz(log_cyc,hist_cyc_log);
cdf_cyc = cdf_cyc/cdf_cyc(end);

% BF distribution
meanBF = mean(stats.BF);
% lognormal
log_BF = log(stats.BF);
logBF = fitdist_trunc(log_BF,@normpdf,@normcdf,log(results.bump_f),[mean(log_BF),std(log_BF)]);
% histogram
rng(0);
burst_freq_ctn = stats.BF+(rand(size(stats.BF))-0.5)*stats.df;   % continuous value
[hist_freq,frq] = hist(burst_freq_ctn,nbins);
[hist_freq_log,log_freq] = hist(log(burst_freq_ctn),nbins);
idx = hist_freq_log==0;	hist_freq_log(idx) = 1;
hist_freq_log = nvalid/(nvalid+sum(idx))*hist_freq_log;
% empirical cdf of BF
df = binwidth(log_freq);
freq_edge = logBF(1)+4*logBF(2)*[-1,1];
log_frq_l = fliplr(log_freq(1)-df:-df:freq_edge(1));
log_frq_r = log_freq(end)+df:df:freq_edge(2);
log_frq = [log_frq_l,log_freq,log_frq_r];
pdf_frq_empr = hist_freq_log/(nvalid*df)*diff(normcdf(log_freq([1,end])+df/2*[-1,1],logBF(1),logBF(2)));
pdf_frq = [normpdf(log_frq_l,logBF(1),logBF(2)),pdf_frq_empr,normpdf(log_frq_r,logBF(1),logBF(2))];
cdf_frq = cumtrapz(log_frq,pdf_frq);
cdf_frq = cdf_frq/cdf_frq(end);

% caculate covariance
X = [stats.AP,stats.CN,stats.BF];
CorrLin = corrcoef(X);
CorrLog = corrcoef(log(X));

% Add results
varnames = {'meanAP','gamAP','expAP','logAP','meanCN','logCN', ...
    'log_cyc','cdf_cyc','meanBF','logBF','log_frq','cdf_frq', ...
    'pk_amp','hist_pks','rd_pks','pdf_pks','CorrLin','CorrLog'};
for i = 1:numel(varnames),	results.(varnames{i})=eval(varnames{i});	end

% Display results
if p.Results.display
    if ~isempty(p.Results.msgfcn)
        dispfcn = p.Results.msgfcn;
    else
        dispfcn = @disp;
    end
    dispfcn(num2str([meanAP,threshold,mean(stats.AP)],...
        ['Mean AS-AP: %.2f ',char(956),'V, AS amplitude threshold: %.2f ',...
        char(956),'V, Mean valid burst peak amplitude: %.2f ',char(956),'V']));
    dispfcn([num2str(expAP,['Exponential mean: %.3f ',char(956),'V']),...
        num2str(logAP,', Lognormal mean/(std): %.3f/(%.3f)'),...
        num2str(gamAP,[', Gamma (k,',char(952),'): (%.3f,%.3f)'])]);
    dispfcn(num2str([meanCN,logCN],...
        'Mean cycle number: %.2f, Lognormal mean/(std): %.3f/(%.3f)'));
    dispfcn(num2str([meanBF,logBF],...
        'Mean burst frequency: %.2f Hz, Lognormal mean/(std): %.3f/(%.3f)'));
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
    if with_lo
        nbins_lo = floor(2*nbins*threshold_dist/maxpeak);
        with_lo = nbins_lo>=2;
    end
    % linear scale
    if with_lo
        i_hi_peaks = peaks>=threshold_dist;
        logAP_hi = fitdist_trunc(log_AP,@normpdf,@normcdf,[log(threshold_dist),Inf], ...
            [mean(log_AP(i_hi_peaks)),std(log_AP(i_hi_peaks))]);
        npk = sum(i_hi_peaks);
        d_pk_amp = threshold_dist/nbins_lo;
        pk_amp_lo = linspace(0,threshold_dist,nbins_lo+1);
        pk_amp_lo = (pk_amp_lo(1:end-1)+pk_amp_lo(2:end))/2;
        pk_amp_hi = pk_amp_lo(end)+d_pk_amp:d_pk_amp:maxpeak;
        h_pks_lo = hist(peaks(~i_hi_peaks),pk_amp_lo);
        h_pks_hi = hist(peaks(i_hi_peaks),pk_amp_hi);
        pk_amp = [pk_amp_lo,pk_amp_hi];
        pk_pdf_scl = npeak*d_pk_amp;
    end
    threshold_idx = find(pk_amp>=threshold,1);
    % log scale
    if plot_logscale
        log_pks_hi = linspace(log(threshold_dist),log(maxpeak),nbins+1);
        log_pks_hi = (log_pks_hi(1:end-1)+log_pks_hi(2:end))/2;
        d_pk_log = binwidth(log_pks_hi);
        log_pks_lo = fliplr(log_pks_hi(1)-d_pk_log:-d_pk_log:log(min(peaks)));
        log_pks = [log_pks_lo,log_pks_hi];  log_pks_exp = exp(log_pks);
        if with_lo
            h_pks_log_hi = hist(log_AP(i_hi_peaks),log_pks_hi);
            h_pks_log_lo = hist(log_AP(~i_hi_peaks),log_pks_lo);
        else
            h_pks_log = hist(log_AP,log_pks);
        end
        threshold_log_idx = find(log_pks>=log(threshold),1);
    end
    if axflag
        figure;
        if plot_logscale,	subplot(211);   end
        ax(1) = gca;
    end
    hold(ax(1),'on');
    if with_lo
        bar(ax(1),pk_amp_hi,h_pks_hi/pk_pdf_scl,'hist');
        set(bar(ax(1),pk_amp_lo,h_pks_lo/pk_pdf_scl,'hist'),'FaceColor','b');
        sig_lo = fitdist_trunc(peaks,@maxwboltzpdf,@maxwboltzcdf, ...
            [0,threshold_dist],(pi/8)^0.5*mean(peaks(~i_hi_peaks)));
        plot(ax(1),pk_amp,(1-npk/npeak)/maxwboltzcdf(threshold_dist,sig_lo)*maxwboltzpdf(pk_amp,sig_lo),'c');
        plot(ax(1),pk_amp(threshold_idx:end),stats.ndetect/npeak*...
            exppdf(pk_amp(threshold_idx:end),expAP)/(1-expcdf(threshold,expAP)),'m','LineWidth',2);
        plot(ax(1),pk_amp,npk/npeak/(1-logncdf(threshold_dist,logAP_hi(1),logAP_hi(2)))...
            *lognpdf(pk_amp,logAP_hi(1),logAP_hi(2)),'y');
    else
        bar(ax(1),pk_amp,pdf_pks,'hist');
    end
    plot(ax(1),pk_amp,gampdf(pk_amp,gamAP(1),gamAP(2)),'r');
    xlim(ax(1),[0,maxpeak]);	xlabel(ax(1),'AS-AP amplitude (\muV)');
    ylabel(ax(1),'Probability/\muV');
    if with_lo
        leg = legend(ax(1),{'High amplitude','Low amplitude',...
            'Maxwell-Boltzmann distribution','Exponential distribution',...
            'Lognormal distribution','Gamma distribution'},'Location','NorthEast');
    else
        leg = legend(ax(1),{'Actual histogram','Gamma distribution'},'Location','NorthEast');
    end
    set(leg,'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(4) = gca;	end
        hold(ax(4),'on');
        if with_lo
            set(bar(ax(4),log_pks_lo,h_pks_log_lo/(npeak*d_pk_log),'hist'),'FaceColor','b');
            bar(ax(4),log_pks_hi,h_pks_log_hi/(npeak*d_pk_log),'hist');
            plot(ax(4),log_pks,(1-npk/npeak)/maxwboltzcdf(threshold_dist,sig_lo)...
                *log_pks_exp.*maxwboltzpdf(log_pks_exp,sig_lo),'c');
            plot(ax(4),log_pks,npk/npeak/(1-logncdf(threshold_dist,logAP_hi(1),logAP_hi(2)))*...
                normpdf(log_pks,logAP_hi(1),logAP_hi(2)),'y');
            plot(ax(4),log_pks(threshold_log_idx:end),stats.ndetect/npeak*log_pks_exp(threshold_log_idx:end)...
                .*exppdf(log_pks_exp(threshold_log_idx:end),expAP)/(1-expcdf(threshold,expAP)),'m','LineWidth',2);
        else
            bar(ax(4),log_pks,h_pks_log/(npeak*d_pk_log),'hist');
        end
        plot(ax(4),log_pks,log_pks_exp.*gampdf(log_pks_exp,gamAP(1),gamAP(2)),'r');
        xlim(ax(4),log_pks([1,end]));	xlabel(ax(4),'AS-AP amplitude (log scale)');
        ylabel(ax(4),'Probability density');
    end
    % CN
    if axflag
        figure;
        if plot_logscale,	subplot(211);   end
        ax(2) = gca;
    end
    hold(ax(2),'on');
    bar(ax(2),cyc,hist_cyc/(nvalid*binwidth(cyc)),'hist');
    plot(ax(2),cyc,lognpdf(cyc,logCN(1),logCN(2)),'r');
    xlabel(ax(2),'Cycle number');	ylabel(ax(2),'Probability density');
    set(legend(ax(2),{'Actual histogram','Lognormal distribution'}),'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(5) = gca;	end
        hold(ax(5),'on');
        bar(ax(5),log_cyc,hist_cyc_log/(nvalid*binwidth(log_cyc)),'hist');
        plot(ax(5),log_cyc,normpdf(log_cyc,logCN(1),logCN(2)),'r');
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
    bar(ax(3),frq,hist_freq/(nvalid*binwidth(frq)),'hist');
    plot(ax(3),frq,lognpdf(frq,logBF(1),logBF(2)),'r');
    xlabel(ax(3),'Burst frequency (Hz)');	ylabel(ax(3),'Probability/Hz');
    set(legend(ax(3),{'Actual histogram','Lognormal distribution'}),'color','none');
    if plot_logscale
        if axflag,	subplot(212);	ax(6) = gca;	end
        hold(ax(6),'on');
        bar(ax(6),log_freq,hist_freq_log/(nvalid*binwidth(log_freq)),'hist');
        plot(ax(6),log_freq,normpdf(log_freq,logBF(1),logBF(2)),'r');
        h = plot(ax(6),log_frq,pdf_frq,'c');
        xlabel(ax(6),'Burst frequency (log scale)');	ylabel(ax(6),'Probability density');
        set(legend(ax(6),h,'Empirical fitting'),'color','none');
    end
end

end
