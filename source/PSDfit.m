function [a,b,outliers,fit_range] = PSDfit( X,Y,freq_band,varargin )
p = inputParser;
validVector = @(x) validateattributes(x,{'numeric'},{'vector'});
validRangePair = @(x) validateattributes(x,{'numeric'},{'nonnegative','integer','numel',2});
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'positive','scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'X',validVector); % frequency
addRequired(p,'Y',validVector); % power density
addRequired(p,'freq_band',validRangePair); % index of signal frequency band
addOptional(p,'tDB',1,validPositiveScalar);	% decibel threshold
addParameter(p,'autofit',false,validTF);	% find fit range automatically
addParameter(p,'tightmode',false,validTF);	% further shrink fit range within freq_band with autofit
addParameter(p,'plot',false,validTF);	% plot iterations
addParameter(p,'framerate',8,validPositiveScalar);	% frame rate for plot
parse(p,X,Y,freq_band,varargin{:})

% inputs
autofit = p.Results.autofit;
plt = p.Results.plot;
freq_band = p.Results.freq_band(:)';
tDB = p.Results.tDB;
tDB = log(db2pow(tDB));
X = log(X(:));
Y = log(Y(:));

% outputs
ind = 1:numel(X);
outliers = [];
P = [0,0];

% fit
if plt
    fig = figure;   hold on;    hh = zeros(1,2);
    hh(1) = plot(X,Y,'b.');
    ylim('auto');	yl = get(gca,'YLim');
    hh(2) = plot(X(freq_band(1))*[1,1],yl,'--','color','#EDB120');
    plot(X(freq_band(2))*[1,1],yl,'--','color','#EDB120');
    axis tight; ylim(yl);
    xlabel('logscale (Hz)');    ylabel('logscale PSD');
end
if autofit
    if plt
        h = -ones(1,4);
        dt = 1/p.Results.framerate;
        it = 0;
    end
    tightmode = p.Results.tightmode;
    bound = freq_band;
    err = [];
    flag = true;
    while flag
        x = X(ind); y = Y(ind);
        FIT;
        if plt
        P0 = P;
        outliers0 = outliers;
        end
        if tightmode
            segs = ind2seg(outliers);   % target outliers
            if ~isempty(segs)
                bound = ind(segs([1,end]));
            end
        end
        outliers_u = err>tDB&~outliers; % upward undesired outliers
        outliers_d = err<-tDB;  % downward outliers
        segs = [ind2seg(outliers_u);ind2seg(outliers_d)];   % undesired outliers
        flag = ~isempty(segs);
        if flag
            segs = [min(segs,[],1);max(segs,[],1)];
            i_l = 1;	i_r = numel(ind);
            c = zeros(1,2); % conditions on each side: 0-not exist, 1-not on edge; 2-on edge
            if ind(segs(1,2))<bound(1)
                c(1) = 1+(i_l==segs(1,1));
            end
            if ind(segs(2,1))>bound(2)
                c(2) = 1+(i_r==segs(2,2));
            end
            flag = any(c);
            if flag
                if c(1)==c(2)
                    % prioritize the one with greater error
                    s = max(abs(err(segs(1,1):segs(1,2))))<max(abs(err(segs(2,1):segs(2,2))));
                    % prioritize the further one
%                     s = X(bound(1))-x(segs(1,2))<x(segs(end,1))-X(bound(2));
                else
                    s = c(1)<c(2);
                end
                switch c(1+s)+2*s
                    case 1
                        i_l = ceil((i_l+segs(1,1))/2);
                    case 2
                        i_l = segs(1,2)+1;
                    case 3
                        i_r = floor((i_r+segs(2,2))/2);
                    case 4
                        i_r = segs(2,1)-1;
                end
                ind = ind(i_l:i_r);
            end
        end
        if ~flag
            FIT;
        end

        if plt
            it = it+1;
            PLOT;
            title(['Iteration # ',num2str(it)]);
            pause(dt);
        end
    end
else
    x = X;	y = Y;
    flag = false;
    FIT;
    if plt
        PLOT;
    end
end

% results
outliers = ind(ind2seg(outliers));
a = exp(P(2));
b = P(1);
fit_range = ind([1,end]);

    function FIT()
        if flag
            maxit = 2;
        else
            maxit = numel(x);
        end
        outliers_pre = false(size(x));
        for k = 1:maxit
            P = polyfit(x(~outliers_pre),y(~outliers_pre),1);
            err = y-(P(1)*x+P(2));
            outliers = err>tDB;
            segs0 = ind2seg(outliers);
            for j = find(ind(segs0(:,2))<freq_band(1) | ind(segs0(:,1))>freq_band(2))
                outliers(segs0(j,1):segs0(j,2)) = false;
            end
            if all(outliers==outliers_pre)
                break;
            else
                outliers_pre = outliers;
            end
        end
    end

    function PLOT()
        figure(fig);
        H = hh;
        labels = {'smoothed PSD','signal frequency band'};
        if autofit
            for j = 1:4
                if isgraphics(h(j))
                    delete(h(j));
                end
            end
            h(1) = plot(x,P0(1)*x+P0(2),'m');
            h(2) = plot(x,y,'c.');
            outlier_out = outliers_u|outliers_d;
            if any(outlier_out)
                h(3) = plot(x(outlier_out),y(outlier_out),'m.');
            else
                h(3) = plot(nan,nan,'m.');
            end
            if any(outliers0)
                h(4) = plot(x(outliers0),y(outliers0),'r.');
            else
                h(4) = plot(nan,nan,'r.');
            end
            H = [H,h];
            labels = [labels,'second fit','normal points' ...
                'undesired outliers','target outliers'];
        end
        if ~flag && any(outliers)
            hhh = zeros(1,2);
            hhh(1) = plot(x,P(1)*x+P(2),'g');
            hhh(2) = plot(x(outliers),y(outliers),'g.');
            H = [H,hhh];
            labels = [labels,'final fit','final target outliers'];
        end
        legend(H,labels,'Location','southwest');
    end

end
