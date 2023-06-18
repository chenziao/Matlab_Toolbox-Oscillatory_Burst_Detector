classdef SynthParam < handle
    % Synthesizing Parameters Class
    properties (Dependent)
        fs                  % synthetic signal sampling rate (Hz)
        dt                  % time step (sec). default: 0.001
        syn_len             % synthetic signal length (sec)
        amp_dist_type;
    end
    properties (SetAccess=private, Hidden)
        Fs = 1000;
        Dt = 1/1000;
        Syn_len = 1000;
        NT                  % total number of time steps
        mu                  % mean of log gaussian copula
        sigma               % std of log gaussian copula
        sigcov              % covariance matrix
        fig_AP_dist         % figure for peak distribution
    end
    properties
        Amp_Dist_Type = 'gamma';    % type of amplitude distribution. {'gamma'},'exponential','lognormal'
        bl_power_match = false; % match only band limited power. default: 0
        power_match = false;	% method for matching the burst power. {0:estimated}, 1:exact
        correlated = true;	% generate correlated attributes AP, CN, BF. default: 1
        empr_CN = true;     % use emperical distribution of cycle number. default: 1
        empr_BF = true;     % use emperical distribution of burst frequency. default: 1
        width_edge = 3;     % burst atom width, multiple of sigma of gaussian envelope. default: 3
        dur_edge            % sigma of gaussian corresponding to stop_perc
        unit_power          % energy of each burst atom with unit amplitude/1 sigma duration
        burst_energy        % total energy of burst trace
        butter_order = [];	% butterworth filter order. default: 6
        afilt               % butterworth filter coefficient vector of a
        bfilt               % butterworth filter coefficient vector of b
        sos                 % butterworth filter second-order section
        g                   % butterworth filter overall system gain
        eval_param          % parameters for amplitude distribution evaluation
        mem_presever = 1e9;	% memory preserved after allocating burst atoms cache
        maxnumsampperdim = 128;	% maximum number of samples to generate along each dimension of burst atom
        minnumburst = 1000; % minimum number of burst atoms to be generated
        randseed            % random seed for synthesizing
        filepath            % file path for the real data file
        rdat_raw            % raw real data information
        rdat                % real data information (adjusted)
        cache               % cache for intermediate synthesized data
        amp_param           % amplitude distribution parameters
        bg_pow_scale        % background power scaling factor
        div_value           % amplitude optimization divergence value
        dist_snapshot       % amplitude distribution snapshot after optimization
        clr_type            % color of each amplitude distribution type for plot
        ks_test_pvalue      % p-value of Kolmogorov-Smirnov test
        bg_pow_scale_whole = false; % apply to whole frequency range when scaling background power. default: 0
    end
    
    methods
        %% Initialize
        function obj = SynthParam(filepath)
            obj.syn_len = obj.Syn_len;
            obj.eval_param = struct();
            obj.fig_AP_dist = struct('fig',[],'h',[0,0],'ylim',[]);
            obj.amp_param = struct('gamma',[],'exponential',[],'lognormal',[]);
            obj.div_value = obj.amp_param;
            obj.bg_pow_scale = struct('gamma',1,'exponential',1,'lognormal',1);
            obj.dist_snapshot = obj.amp_param;
            obj.clr_type = struct('gamma','b','exponential','y','lognormal','g');
            obj.ks_test_pvalue = obj.amp_param;
            if nargin>0
                obj.filepath = filepath;
                setup(obj);
            end
        end
        %% Get/Set dependent properties
        function val = get.fs(obj)
            val = obj.Fs;
        end
        function set.fs(obj,val)
            obj.Fs = val;
            obj.Dt = 1/val;
            obj.NT = ceil(obj.Syn_len/obj.Dt);
        end
        function val = get.dt(obj)
            val = obj.Dt;
        end
        function set.dt(obj,val)
            obj.Dt = val;
            obj.Fs = 1/val;
            obj.NT = ceil(obj.Syn_len/obj.Dt);
        end
        function val = get.syn_len(obj)
            val = obj.Syn_len;
        end
        function set.syn_len(obj,val)
            obj.Syn_len = val;
            obj.NT = ceil(obj.Syn_len/obj.Dt);
        end
        function val = get.amp_dist_type(obj)
            val = obj.Amp_Dist_Type;
        end
        function set.amp_dist_type(obj,val)
            if ~isfield(obj.amp_param,val)
                warndlg(['Amplitude distribution type not defined. ', '''', ...
                    obj.Amp_Dist_Type,''' distribution will be used.'],'Value Warning');
                return;
            end
            obj.Amp_Dist_Type = val;
            if isempty(obj.rdat),   obj.loadrealdata;   end
            obj.scale_bg_pow;
        end
        
        %% Setup
        % Load real data from given file.
        function obj = loadrealdata(obj)
            if isempty(obj.filepath)
                error('Property "filepath" needs to be specified.');
            end
            obj.rdat = load(obj.filepath,'-mat');
            obj.rdat_raw = obj.rdat;
            obj.scale_bg_pow;
        end
        % Setup Synthesization parameters
        function [obj,flag] = setup(obj,optimize)
            % Load real data from given file.
            if isempty(obj.rdat),   obj.loadrealdata;   end
            % Select raw real data for optimization
            if nargin<2 || isempty(optimize) || ~optimize
                RD = obj.rdat;
            else
                RD = obj.rdat_raw;
            end
            % Initialize figure
            if ~isempty(obj.fig_AP_dist.fig) && ishandle(obj.fig_AP_dist.fig)
                obj.fig_AP_dist = struct('fig',[],'h',[0,0],'ylim',[]);
            end
            % Setup synthesization parameters.
            if obj.Fs>RD.fs
                warndlg(['Using sampling rate greater than the in vivo data '...
                    'will cause extrapolating the power spectral density '...
                    'beyond the Nyquist frequency of the in vivo data.'],'Value Warning');
            end
            if obj.Fs<=RD.sig_range_f(2)*2
                errordlg(num2str(RD.sig_range_f(2)*2,['Sampling rate must be ' ...
                    'greater than twice the upper bound of the signal frequency ' ...
                    'range: %g Hz.']),'Value Error');
                flag = 0;   return;
            end
            obj.mu = [RD.logn_AP(1),RD.logn_CN(1),RD.logn_BF(1)];
            obj.sigma = [RD.logn_AP(2),RD.logn_CN(2),RD.logn_BF(2)];
            obj.sigcov = RD.CorrLog;
            obj.dur_edge = fminbnd(@(x) abs(normpdf(x)-RD.stop_perc*normpdf(0)),0,3);
            f = RD.f;
            PSD_smoo = RD.PSD_smoo;
            Fs2 = obj.Fs/2;
            if obj.bl_power_match
                bg_range_i = RD.sig_range_i(1):RD.sig_range_i(2);
            else
                bg_range_i = RD.bg_range_i;
                if Fs2<f(bg_range_i(end))
                    warndlg(['Sampling rate too low. Ignoring part of the ', ...
                        'signal power beyond the Nyquist frequency.'],'Value Warning');
                    bg_range_i = [bg_range_i(f(bg_range_i)<Fs2),numel(f)+1];
                    PSD_smoo(end+1) = interp1(f,PSD_smoo,Fs2);
                    f(end+1) = Fs2;
                end
            end
            obj.burst_energy = trapz(f(bg_range_i),PSD_smoo(bg_range_i)- ...
                getBackgroundFit(RD,f(bg_range_i)))*obj.Syn_len;
            if obj.burst_energy <=0
                errordlg(['Total power is lower than the background power, ', ...
                    'resulting in negative signal power.'],'Characterization Error');
                flag = 0;   return;
            end
            obj.unit_power = sqrt(pi)*erf(obj.width_edge)/4;
            if isempty(obj.butter_order)
                obj.butter_order = RD.butter_order;
            end
%             [obj.bfilt,obj.afilt] = butter(obj.butter_order,RD.sig_range_f*2/obj.Fs);
            [Z,P,K] = butter(obj.butter_order,RD.sig_range_f*2/obj.Fs);
            [obj.sos,obj.g] = zp2sos(Z,P,K);	% using sos for numerical stability
            [obj.bfilt,obj.afilt] = zp2tf(Z,P,K);
            flag = 1;
        end
        % Setup amplitude distribution evaluation parameters
        function obj = setup_eval_amp_dist(obj,div_type,fit_AP,alpha,threshold,randburst,verbose,plotdist)
            % div_type is Kullback-Leibler divergence or Jensen-Shannon divergence. {'KL'},'JS'.
            % fit_AP is a flag to use amplitude peak as target for fitting. defualt: false, using amplitude.
            % alpha is a value between 0 and 1 added to the count for each event,
            % so called add-alpha smoothing, also called Laplace smoothing. default:1
            % threshold is a flag for considering only amplitude higher than a
            % threshold when evaluating distribution divergence. default: false
            % randburst is a flag for using random permuted burst. default: true
            % verbose is the flag for displaying details such as progress. default: false
            % plotdist is the flag for plotting distribution after each evaluation. default: false
            RD = obj.rdat_raw;
            if nargin<2 || ~ischar(div_type)
                obj.eval_param.div_type = 'KL';
            else
                obj.eval_param.div_type = div_type;
            end
            if nargin<3 || ~isscalar(fit_AP)
                obj.eval_param.fit_AP = false;
            else
                obj.eval_param.fit_AP = logical(fit_AP);
            end
            if nargin<4 || ~isscalar(alpha)
                alpha = 1;
            end
            if nargin<5 || ~isscalar(threshold) || ~threshold
                threshold_idx = 1;
            else
                threshold_idx = find(RD.AP_edges>=RD.threshold,1);
            end
            if obj.eval_param.fit_AP
                desired_hist = RD.hist_AP(threshold_idx:end);
            else
                desired_hist = RD.hist_amp(threshold_idx:end);
            end
            switch div_type
                case 'KL'
                    obj.eval_param.divfun = @(x) KLDiv(desired_hist,x,alpha);
                case 'JS'
                    obj.eval_param.divfun = @(x) JSDiv(desired_hist,x,alpha);
                otherwise
                    obj.eval_param.divfun = [];
                    error('"div_type" must be "KL" or "JS".');
            end
            obj.eval_param.threshold_idx = threshold_idx;
            if nargin<6 || ~isscalar(randburst)
                obj.eval_param.randburst = true;
            else
                obj.eval_param.randburst = logical(randburst);
            end
            if nargin<7 || ~isscalar(verbose)
                obj.eval_param.verbose = false;
            else
                obj.eval_param.verbose = logical(verbose);
            end
            if nargin<8 || ~isscalar(plotdist)
                obj.eval_param.plotdist = false;
            else
                obj.eval_param.plotdist = logical(plotdist);
            end
        end
        % Reset to initial random state before evaluation
        function obj = reset_randstate(obj)
            if ~isfield(obj.cache,'initrandstate')
                error('Cache for intermediate synthesized data not created yet.');
            else
                obj.cache.randstate = obj.cache.initrandstate;
            end
        end
        % Scale background power after optimization
        function obj = scale_bg_pow(obj)
            s = obj.bg_pow_scale.(obj.Amp_Dist_Type);
            obj.rdat.tDB = obj.rdat_raw.tDB-pow2db(s);
            obj.rdat.fit_a = s*obj.rdat_raw.fit_a;
            if isfield(obj.rdat,'aperiodic_params')
                obj.rdat.aperiodic_params(1) = obj.rdat.aperiodic_params(1)+log10(s);
            end
            obj.rdat.pow_bg = s*obj.rdat_raw.pow_bg;
            obj.rdat.SNR = (obj.rdat_raw.SNR+1)/s-1;
            obj.rdat.pow_sig = obj.rdat.pow_bg*obj.rdat.SNR;
            obj.rdat.pow_bg_filt = s*obj.rdat_raw.pow_bg_filt;
            obj.rdat.SNR_filt = (obj.rdat_raw.SNR_filt+1)/s-1;
            obj.rdat.pow_sig_filt = obj.rdat.pow_bg_filt*obj.rdat.SNR_filt;
        end
        % Save object to file
        function obj = saveobj(obj)
            obj.cache = [];    obj.rdat = [];
            obj.fig_AP_dist = struct('fig',[],'h',[0,0],'ylim',[]);
        end
        
        %% Functions
        % Calculate background PSD.
        function [PSD_bg,f,Nfft,bg_range_i,PSD_tot] = psd_background(obj,Fs)
            if nargin<2,    Fs = obj.Fs;	end
            Nfft = 2*round(Fs/2/binwidth(obj.rdat.f));
            f = linspace(0,Fs/2,Nfft/2+1);
            f_len = find(f<=obj.rdat.f(end),1,'last');
            PSD_fit = getBackgroundFit(obj.rdat,f);
            PSD_tot = [interp1(obj.rdat.f,obj.rdat.PSD_smoo, ...
                f(1:f_len)),PSD_fit(f_len+1:Nfft/2+1)];
            bg_range_i = [find(f<=obj.rdat.bg_range_f(1),1,'last'), ...
                find(f>=obj.rdat.bg_range_f(2),1,'first')];
            if length(bg_range_i)==1,	bg_range_i(2) = Nfft/2+1;	end
            bg_range_i = bg_range_i(1):bg_range_i(2);
            if obj.bg_pow_scale_whole
                PSD_tot_i = PSD_tot(bg_range_i);
                PSD_tot = obj.bg_pow_scale.(obj.Amp_Dist_Type)*PSD_tot;
                PSD_tot(bg_range_i) = PSD_tot_i;
            end
            PSD_bg = PSD_tot;
            PSD_bg(bg_range_i) = PSD_fit(bg_range_i);
            PSD_bg(1) = 0;
        end
        % Generate background trace
        function backgroundtrace = gen_background(obj,dist_type)
            if nargin<2,	dist_type = 'Gaussian';	end
            if ~isempty(obj.randseed),	rng(obj.randseed+1);	end	% Set random seed
            [PSD_bg,f,~,bg_range_i] = obj.psd_background;
            filt_n = min(floor(f(end)),1000);   % 1 Hz resolution
            backgroundtrace = noiseGenPSD(obj.NT,PSD_bg,f,filt_n,dist_type,0,f(bg_range_i([1,end])));
            obj.cache.backgroundTrace = backgroundtrace;
        end
        % Calculate numerical PDF of gaussian copula.
        function [Xg,PX] = copula_pdf(obj,nsig,dsig)
            % Xg is the meshgrid of the 3 properties of burst atoms
            % PX is the probability of each burst atom with properties Xg
            % nsig is the cutoff range of normal distribution to include.
            % dsig is the width of each grid in the meshgrid.
            if nargin<2||~isscalar(nsig),	nsig = 4.5;	end
            if nargin<3||~isscalar(dsig),	dsig = 0.2;	end
            gv = -nsig:dsig:nsig;
            ng = length(gv)^3;  Xg = zeros(ng,4);
            [Xg(ng+1:2*ng),Xg(1:ng),Xg(2*ng+1:3*ng)] = meshgrid(gv);
            if obj.correlated
                PX = mvnpdf(Xg(:,1:3),[0,0,0],obj.sigcov);
            else
                PX = zeros(ng,3);
                [PX(ng+1:2*ng),PX(1:ng),PX(2*ng+1:3*ng)] = meshgrid(normpdf(gv));
                PX = prod(PX,2);
            end
            PX = PX/sum(PX);
            if obj.empr_CN
                Xg(:,2) = interp1(obj.rdat.cdf_log_CN,obj.rdat.log_CNs,normcdf(Xg(:,2)));
            else
                Xg(:,2) = Xg(:,2)*obj.sigma(2)+obj.mu(2);
            end
            if obj.empr_BF
                Xg(:,3) = interp1(obj.rdat.cdf_log_BF,obj.rdat.log_BFs,normcdf(Xg(:,3)));
            else
                Xg(:,3) = Xg(:,3)*obj.sigma(3)+obj.mu(3);
            end
            Xg(:,4) = exp(Xg(:,2)-Xg(:,3)); % 4th column stores the duration
            if ~nargout
                obj.cache.Xg = Xg;
                obj.cache.PX = PX;
            end
        end
        % Calculate proportion of band limited power of busrt atom
        function p = bl_power_prop(obj,bf,dur)
            % p is the proportion of power of a burst atom within the band
            % bf is the burst frequency (column vector)
            % dur is the duration between +/-1 sigma (column vector)
            sig_range_f = obj.rdat.sig_range_f;
            bw = (2^-.5/pi)./dur;
            p = normcdf(sig_range_f(2),bf,bw)-normcdf(sig_range_f(1),bf,bw);
        end
        % Plot amplitude peak distribution after each evaluation
        function [h1,h2] = plot_AP_dist(obj,hist_AP,ax,h2)
            if nargin<3
                newfig = isempty(obj.fig_AP_dist.fig) || ~ishandle(obj.fig_AP_dist.fig);
                if newfig
                    obj.fig_AP_dist.fig = figure;    hold on;
                    xlabel(['Amplitude peak (',obj.rdat.unit,')']);
                    ylabel(['Rate density (Hz/',obj.rdat.unit,')']);
                    xlim(obj.rdat.AP_edges([1,end]));
                else
                    figure(obj.fig_AP_dist.fig);	hold on;
                end
                if isequal(obj.fig_AP_dist.h(1),0) || ~ishandle(obj.fig_AP_dist.h(1))
                    obj.fig_AP_dist.h(1) = plot(obj.rdat.AP_ctrs,obj.rdat.rd_AP,'r');
                    obj.fig_AP_dist.ylim = get(gca,'ylim');
                end
                if ~isequal(obj.fig_AP_dist.h(2),0) && ishandle(obj.fig_AP_dist.h(2))
                    delete(obj.fig_AP_dist.h(2));
                end
                obj.fig_AP_dist.h(2) = plot(obj.rdat.AP_ctrs,...
                    hist_AP/obj.Syn_len/binwidth(obj.rdat.AP_ctrs),'b');
                ylim(obj.fig_AP_dist.ylim);
                if newfig
                    legend(obj.fig_AP_dist.h,{'Observed','Synthetic'},'Location','NorthEast');
                end
                h1 = [];	h2 = [];
            else
                hold(ax,'on');
                if isempty(obj.fig_AP_dist.fig) || ~ishandle(obj.fig_AP_dist.fig)
                    obj.fig_AP_dist.fig = ax;	
                    xlabel(ax,['Amplitude peak (',obj.rdat.unit,')']);
                    ylabel(ax,['Rate density (Hz/',obj.rdat.unit,')']);
                    xlim(ax,obj.rdat.AP_edges([1,end]));
                end
                if isequal(obj.fig_AP_dist.h(1),0) || ~ishandle(obj.fig_AP_dist.h(1))
                    obj.fig_AP_dist.h(1) = plot(ax,obj.rdat.AP_ctrs, ...
                    obj.rdat.rd_AP,'color',0.5*[1,1,1],'Linewidth',2);
                    obj.fig_AP_dist.ylim = get(ax,'ylim');
                end
                if ~isequal(h2,0) && ishandle(h2)
                    delete(h2);
                end
                h1 = obj.fig_AP_dist.h(1);
                h2 = plot(ax,obj.rdat.AP_ctrs,hist_AP/obj.Syn_len/ ...
                    binwidth(obj.rdat.AP_ctrs),':','Linewidth',2);
                ylim(ax,obj.fig_AP_dist.ylim);
            end
        end
        % Record divergence value after optimization
        function obj = record_div_val(obj,div_val)
            obj.div_value.(obj.amp_dist_type) = div_val;
        end
        % Plot amplitude and amplitude peak distribution of saved evaluation
        function [fig,h1,h2] = plot_eval_dist(obj,amp_dist_type,single_plot,ax,h1,h2)
            if obj.eval_param.fit_AP
                fit_target = 'amplitude peak';
            else
                fit_target = 'amplitude';
            end
            if nargin<3 || isempty(single_plot)
                single_plot = false;
            end
            if nargin<4
                fig = figure;   ax = gca;
                clr = 'r';  LS = '-';   LW = 1;
            else
                fig = [];
                clr = 0.5*[1,1,1];  LS = ':';   LW = 2;
                single_plot = true;
                if ~isequal(h1,0) && ishandle(h1)
                    delete(h1);
                end
                if ~isequal(h2,0) && ishandle(h2)
                    delete(h2);
                end
            end
            if ~single_plot || ~obj.eval_param.fit_AP
                if ~single_plot,    subplot(211);   ax = gca;   end
                hold(ax,'on');
                h1 = plot(ax,obj.rdat.amp_ctrs,obj.rdat.pdf_amp,'color',clr,'Linewidth',LW);
                h2 = plot(ax,obj.rdat.amp_ctrs,obj.dist_snapshot.(amp_dist_type){1}, ...
                    LS,'color',obj.clr_type.(amp_dist_type),'Linewidth',LW);
                xlabel(ax,['Amplitude (',obj.rdat.unit,')']);
                ylabel(ax,'Probability Density');
            end
            if ~single_plot || obj.eval_param.fit_AP
                if ~single_plot,    subplot(212);	ax = gca;   end
                hold(ax,'on');
                h1 = plot(ax,obj.rdat.AP_ctrs,obj.rdat.rd_AP,'color',clr,'Linewidth',LW);
                h2 = plot(ax,obj.rdat.AP_ctrs,obj.dist_snapshot.(amp_dist_type){2}, ...
                    LS,'color',obj.clr_type.(amp_dist_type),'Linewidth',LW);
                xlabel(ax,['Amplitude peak (',obj.rdat.unit,')']);
                ylabel(ax,['Rate density (Hz/',obj.rdat.unit,')']);
            end
            if ~single_plot,    subplot(211);   end
            if nargin<4
                title(['Fit ',fit_target,' with ',amp_dist_type,' burst atoms']);
                legend({'Observed','Synthetic'},'Location','NorthEast');
            end
        end
        % Find distribution type with minimum divergence value
        function [dist_type,types] = best_amp_dist_type(obj)
            amp_dist_types = reshape(fieldnames(obj.div_value),1,[]);
            div_values = struct2cell(obj.div_value);
            idx = find(~cellfun(@isempty,div_values));
            types = amp_dist_types(idx);
            if isempty(idx)
                dist_type = [];
                warndlg('Amplitude optimization has not been done.');
            else
                [~,I] = min(cell2mat(div_values(idx)));
                dist_type = amp_dist_types{idx(I)};
            end
        end
    end
    
end