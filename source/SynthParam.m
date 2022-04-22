classdef SynthParam < handle
    % Synthesizing Parameters Class
    properties (Dependent)
        fs                  % synthetic signal sampling rate (Hz).
        dt                  % time step (sec). default: 0.001
        syn_len             % synthetic signal length (sec).
    end
    properties (SetAccess=private, Hidden)
        Fs = 1000;
        Dt = 1/1000;
        Syn_len = 1000;
        NT                  % total number of time steps
        mu                  % mean of log gaussian copula
        sigma               % std of log gaussian copula
        sigcov              % covariance matrix
        fig_pkdist          % figure for peak distribution
    end
    properties
        amp_dist_type = 'gamma';	% Type of amplitude distribution. {'gamma'},'exponential','lognormal'
        power_match = false;	% method for matching the burst power. {0:analytical}, 1:emperical
        correlated = true;	% generate correlated amplitude, frequency, cycles. default: 1
        empr_CN = false;	% use emperical distribution of cycle number. default: 0
        empr_BF = false;	% use emperical distribution of burst frequency. default: 0
        width_edge = 3;     % synthetic burst width in sigma of gaussian. default: 3
        dur_edge            % sigma of gaussian corresponding to perc_edge
        unit_power          % power of each burst atom with unit amplitude
        burst_energy        % total energy of burst trace
        butter_order = [];	% butterworth filter order. default: 6
        afilt               % butterworth filter coefficient vector of a
        bfilt               % butterworth filter coefficient vector of b
        sos                 % butterworth filter second-order section
        g                   % butterworth filter overall system gain
        eval_param          % parameters for amplitude distribution evaluation
        mem_presever = 1e9;	% Memory preserved after allocating burst atoms buffer.
        maxnumsampperdim = 128;	% Maximum number of samples to generate along each dimension of burst atom.
        randseed            % random seed for synthesizing.
        filepath            % file path for the real data file.
        rdat                % real data information.
        buffer              % buffer for intermediate synthesized data.
        amp_param           % amplitude distribution parameters.
        div_value           % amplitude optimization divergence value.
        dist_snapshot       % amplitude distribution snapshot after optimization.
    end
    
    methods
        %% Initialize
        function obj = SynthParam(filepath)
            obj.syn_len = obj.Syn_len;
            obj.eval_param = struct();
            obj.fig_pkdist = struct('fig',[],'h',[0,0],'ylim',[]);
            obj.amp_param = struct('exponential',[],'lognormal',[],'gamma',[]);
            obj.div_value = obj.amp_param;
            obj.dist_snapshot = obj.amp_param;
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
        
        %% Setup
        % Load real data from given file.
        function obj = loadrealdata(obj)
            if isempty(obj.filepath)
                error('Property "filepath" needs to be specified.');
            end
            obj.rdat = load(obj.filepath,'-mat');
        end
        % Setup Synthesization parameters
        function [obj,flag] = setup(obj)
            % Load real data from given file.
            if isempty(obj.rdat),   obj.loadrealdata;   end
            % Initialize figure
            if ~isempty(obj.fig_pkdist.fig) && ishandle(obj.fig_pkdist.fig)
                obj.fig_pkdist = struct('fig',[],'h',[0,0],'ylim',[]);
            end
            % Setup synthesization parameters.
            if obj.Fs>obj.rdat.fs
                warndlg(['Using sampling rate greater than the in vivo data '...
                    'will cause extrapolating the power spectral density '...
                    'beyond the Nyquist frequency of the in vivo data.'],'Value Warning');
            end
            if obj.Fs<=obj.rdat.bump_f(2)*2
                errordlg(num2str(obj.rdat.bump_f(2)*2,['Sampling rate must be ' ...
                    'greater than twice the upper bound of the frequency range ' ...
                    'of signal: %g Hz.']),'Value Error');
                flag = 0;   return;
            end
            obj.mu = [obj.rdat.logAP(1),obj.rdat.logCN(1),obj.rdat.logBF(1)];
            obj.sigma = [obj.rdat.logAP(2),obj.rdat.logCN(2),obj.rdat.logBF(2)];
            obj.sigcov = obj.rdat.CorrLog;
            obj.dur_edge = fminbnd(@(x) abs(normpdf(x)-obj.rdat.stop_perc*normpdf(0)),0,3);
            f = obj.rdat.f;
            fit_i = obj.rdat.fit_i;
            PSD_smoo = obj.rdat.PSD_smoo;
            Fs2 = obj.Fs/2;
            if Fs2<f(fit_i(end))
                warndlg(['Sampling rate too low. Ignoring part of the ', ...
                    'signal power beyond the Nyquist frequency.'],'Value Warning');
                fit_i = [fit_i(f(fit_i)<Fs2),numel(f)+1];
                PSD_smoo(end+1) = interp1(f,PSD_smoo,Fs2);
                f(end+1) = Fs2;
            end
            obj.burst_energy = trapz(f(fit_i),PSD_smoo(fit_i)-...
                obj.rdat.fit_a*f(fit_i).^obj.rdat.fit_b)*obj.Syn_len;
            if obj.burst_energy <=0
                errordlg(['Total power within frequency range of interest ', ...
                    'is lower than the background power.'],'Characterization Error');
                flag = 0;   return;
            end
            obj.unit_power = sqrt(pi)*erf(obj.width_edge)/obj.width_edge/4;
            if isempty(obj.butter_order)
                obj.butter_order = obj.rdat.butter_order;
            end
%             [obj.bfilt,obj.afilt] = butter(obj.butter_order,obj.rdat.bump_f*2/obj.Fs);
            [Z,P,K] = butter(obj.butter_order,obj.rdat.bump_f*2/obj.Fs);
            [obj.sos,obj.g] = zp2sos(Z,P,K);	% using sos to avoid numeric error
            [obj.bfilt,obj.afilt] = zp2tf(Z,P,K);
            flag = 1;
        end
        % Setup amplitude distribution evaluation parameters
        function obj = setup_eval_amp_dist(obj,div_type,fit_ASAP,alpha,threshold,randburst,verbose,plotdist)
            % div_type is Kullback-Leibler divergence or Jensen-Shannon divergence. {'KL'},'JS'.
            % fit_ASAP is a flag to use AS-AP as objective for fitting. defualt: false, using AS amplitude.
            % alpha is a value between 0 and 1 added to the count for each event,
            % so called add-alpha smoothing, also called Laplace smoothing. default:1
            % threshold is a flag for considering only amplitude higher than a
            % threshold when evaluating distribution divergence. default: false
            % randburst is a flag for using random permuted burst. default: true
            % verbose is the flag for displaying details such as progress. default: false
            % plotdist is the flag for plotting distribution after each evaluation. default: false
            if nargin<2 || ~ischar(div_type)
                obj.eval_param.div_type = 'KL';
            else
                obj.eval_param.div_type = div_type;
            end
            if nargin<3 || ~isscalar(fit_ASAP)
                obj.eval_param.fit_ASAP = false;
            else
                obj.eval_param.fit_ASAP = logical(fit_ASAP);
            end
            if nargin<4 || ~isscalar(alpha)
                alpha = 1;
            end
            if nargin<5 || ~isscalar(threshold) || ~threshold
                threshold_idx = 1;
            else
                threshold_idx = find(obj.rdat.pk_amp>=obj.rdat.threshold,1);
            end
            if obj.eval_param.fit_ASAP
                desired_hist = obj.rdat.hist_pks(threshold_idx:end);
            else
                desired_hist = obj.rdat.hist_amp(threshold_idx:end);
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
            if ~isfield(obj.buffer,'initrandstate')
                error('Buffer for intermediate synthesized data not created yet.');
            else
                obj.buffer.randstate = obj.buffer.initrandstate;
            end
        end
        % Save object to file
        function obj = saveobj(obj)
            obj.buffer = [];    obj.rdat = [];
            obj.fig_pkdist = struct('fig',[],'h',[0,0],'ylim',[]);
        end
        
        %% Functions
        % Calculate background PSD.
        function [PSD_bg,f,Nfft,fit_i,PSD_tot] = psd_background(obj,Fs)
            if nargin<2,    Fs = obj.Fs;	end
            Nfft = 2*round(Fs/2/binwidth(obj.rdat.f));
            f = linspace(0,Fs/2,Nfft/2+1);
            f_len = find(f<=obj.rdat.f(end),1,'last');
            PSD_fit = obj.rdat.fit_a*f.^obj.rdat.fit_b;
            PSD_tot = [interp1(obj.rdat.f,obj.rdat.PSD_smoo,f(1:f_len)),PSD_fit(f_len+1:Nfft/2+1)];
            fit_i = [find(f<=obj.rdat.fit_f(1),1,'last'),find(f>=obj.rdat.fit_f(2),1,'first')];
            if length(fit_i)==1,	fit_i(2) = Nfft/2+1;	end
            PSD_bg = PSD_tot;
            PSD_bg(fit_i(1):fit_i(2)) = PSD_fit(fit_i(1):fit_i(2));
            PSD_bg(1) = 0;
        end
        % Generate background trace
        function backgroundtrace = gen_background(obj,dist_type)
            if nargin<2,	dist_type = 'Gaussian';	end
            if ~isempty(obj.randseed),	rng(obj.randseed+1);	end	% Set random seed
            [PSD_fit,f,~,fit_i] = obj.psd_background;
            filt_n = min(floor(f(end)),1000);   % 1 Hz resolution
            if nargout
                backgroundtrace = noiseGenPSD(obj.NT,PSD_fit,f,filt_n,dist_type,0,f(fit_i));
            else
                obj.buffer.backgroundTrace = noiseGenPSD(obj.NT,PSD_fit,f,filt_n,dist_type,0,f(fit_i));
            end
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
                Xg(:,2) = interp1(obj.rdat.cdf_cyc,obj.rdat.log_cyc,normcdf(Xg(:,2)));
            else
                Xg(:,2) = Xg(:,2)*obj.sigma(2)+obj.mu(2);
            end
            if obj.empr_BF
                Xg(:,3) = interp1(obj.rdat.cdf_frq,obj.rdat.log_frq,normcdf(Xg(:,3)));
            else
                Xg(:,3) = Xg(:,3)*obj.sigma(3)+obj.mu(3);
            end
            Xg(:,4) = exp(Xg(:,2)-Xg(:,3)); % 4th column stores the duration
            if ~nargout
                obj.buffer.Xg = Xg;
                obj.buffer.PX = PX;
            end
        end
        % Plot peak amplitude distribution after each evaluation
        function [h1,h2] = plot_peak_dist(obj,hist_pks,ax,h2)
            if nargin<3
                newfig = isempty(obj.fig_pkdist.fig) || ~ishandle(obj.fig_pkdist.fig);
                if newfig
                    obj.fig_pkdist.fig = figure;    hold on;
                    xlabel('AS-AP amplitude (\muV)');	ylabel('Rate density (Hz/\muV)');
                    xlim([0,obj.rdat.pk_amp(end)]);
                else
                    figure(obj.fig_pkdist.fig);	hold on;
                end
                if isequal(obj.fig_pkdist.h(1),0) || ~ishandle(obj.fig_pkdist.h(1))
                    obj.fig_pkdist.h(1) = plot(obj.rdat.pk_amp,obj.rdat.rd_pks,'r');
                    obj.fig_pkdist.ylim = get(gca,'ylim');
                end
                if ~isequal(obj.fig_pkdist.h(2),0) && ishandle(obj.fig_pkdist.h(2))
                    delete(obj.fig_pkdist.h(2));
                end
                obj.fig_pkdist.h(2) = plot(obj.rdat.pk_amp,...
                    hist_pks/obj.Syn_len/binwidth(obj.rdat.pk_amp),'b');
                ylim(obj.fig_pkdist.ylim);
                if newfig
                    legend(obj.fig_pkdist.h,{'In Vivo','Synthetic'},'Location','NorthEast');
                end
                h1 = [];	h2 = [];
            else
                hold(ax,'on');
                if isempty(obj.fig_pkdist.fig) || ~ishandle(obj.fig_pkdist.fig)
                    obj.fig_pkdist.fig = ax;	
                    xlabel(ax,'AS-AP amplitude (\muV)');
                    ylabel(ax,'Rate density (Hz/\muV)');
                    xlim(ax,[0,obj.rdat.pk_amp(end)]);
                end
                if isequal(obj.fig_pkdist.h(1),0) || ~ishandle(obj.fig_pkdist.h(1))
                    obj.fig_pkdist.h(1) = plot(ax,obj.rdat.pk_amp, ...
                    obj.rdat.rd_pks,'color',0.5*[1,1,1],'Linewidth',2);
                    obj.fig_pkdist.ylim = get(ax,'ylim');
                end
                if ~isequal(h2,0) && ishandle(h2)
                    delete(h2);
                end
                h1 = obj.fig_pkdist.h(1);
                h2 = plot(ax,obj.rdat.pk_amp,hist_pks/obj.Syn_len/ ...
                    binwidth(obj.rdat.pk_amp),':','Linewidth',2);
                ylim(ax,obj.fig_pkdist.ylim);
            end
        end
        % Record divergence value after optimization
        function obj = record_div_val(obj,div_val)
            obj.div_value.(obj.amp_dist_type) = div_val;
        end
        % Plot AS amplitude and AS-AP distribution of saved evaluation
        function [fig,h1,h2] = plot_eval_dist(obj,amp_dist_type,single_plot,ax,h1,h2)
            if obj.eval_param.fit_ASAP
                fit_target = 'AS-AP';
            else
                fit_target = 'AS amplitude';
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
            if ~single_plot || ~obj.eval_param.fit_ASAP
                if ~single_plot,    subplot(211);   ax = gca;   end
                hold(ax,'on');
                h1 = plot(ax,obj.rdat.as_amp,obj.rdat.pdf_amp,'color',clr,'Linewidth',LW);
                h2 = plot(ax,obj.rdat.as_amp,obj.dist_snapshot.(amp_dist_type){1},LS,'Linewidth',LW);
                xlabel(ax,'AS Amplitude (\muV)');   ylabel(ax,'Probability Density');
            end
            if ~single_plot || obj.eval_param.fit_ASAP
                if ~single_plot,    subplot(212);	ax = gca;   end
                hold(ax,'on');
                h1 = plot(ax,obj.rdat.pk_amp,obj.rdat.rd_pks,'color',clr,'Linewidth',LW);
                h2 = plot(ax,obj.rdat.pk_amp,obj.dist_snapshot.(amp_dist_type){2},LS,'Linewidth',LW);
                xlabel(ax,'AS-AP amplitude (\muV)');    ylabel(ax,'Rate density (Hz/\muV)');
            end
            if ~single_plot,    subplot(211);   end
            if nargin<4
                title(['Fit ',fit_target,' with ',amp_dist_type,' burst atoms']);
                legend({'In Vivo','Synthetic'},'Location','NorthEast');
            end
        end
        % Find distribution type with minimum divergence value
        function dist_type = best_amp_dist_type(obj)
            amp_dist_types = fieldnames(obj.div_value);
            div_values = struct2cell(obj.div_value);
            idx = find(~cellfun(@isempty,div_values));
            if isempty(idx)
                error('Amplitude optimization has not been done.');
                dist_type = [];
            else
                [~,I] = min(cell2mat(div_values(idx)));
                dist_type = amp_dist_types{idx(I)};
            end
        end
    end
    
end