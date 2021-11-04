function [a,b,outliers,fit_range] = PSDfit( X,Y,freq_range,varargin )
p = inputParser;
validVector = @(x) validateattributes(x,{'numeric'},{'vector'});
validRangePair = @(x) validateattributes(x,{'numeric'},{'nonnegative','integer','numel',2});
validPositiveScalar = @(x) validateattributes(x,{'numeric'},{'positive','scalar'});
validTF = @(x) validateattributes(x,{'logical','numeric'},{'scalar'});
addRequired(p,'X',validVector); % frequency
addRequired(p,'Y',validVector); % power density
addRequired(p,'freq_range',validRangePair); % index of frequency range of interest
addOptional(p,'nDB',1,validPositiveScalar);	% decibel threshold
addParamValue(p,'autofit',false,validTF);	% find fit range automatically
addParamValue(p,'plot',false,validTF);	% plot iterations
addParamValue(p,'framerate',5,validPositiveScalar);	% frame rate for plot
parse(p,X,Y,freq_range,varargin{:})

% inputs
freq_range = p.Results.freq_range(:)';
nDB = p.Results.nDB;
tDB = log(db2pow(nDB));
X = log(X(:));
Y = log(Y(:));

% outputs
ind = 1:numel(X);
outliers0 = [];
outliers = [];
P0 = [0,0];
P = [0,0];

% fit
plt = p.Results.plot;
if plt
    fig = figure;   hold on;
    plot(X,Y,'b.');
    ylim('auto');	yl = get(gca,'YLim');
    plot(X(freq_range(1))*[1,1],yl,'y:');
    plot(X(freq_range(2))*[1,1],yl,'y:');
    axis tight; ylim(yl);
    xlabel('logscale (Hz)');    ylabel('logscale PSD');
    h = -ones(1,4);
end
if p.Results.autofit
    if plt
        dt = 1/p.Results.framerate;
        it = 0;
    end
    flag = true;
    bound = [freq_range;freq_range];
    c0 = [false(1,4),true];
    while flag
        x = X(ind); y = Y(ind);
        FIT;
        err = y-(P0(1)*x+P0(2));
        outliers_l = err<-tDB;
        outliers_h = err>tDB;
        segs = ind2seg(outliers_h);
        bound(2,:) = freq_range;
        for i = find(ind(segs(:,1))<=freq_range(2) & ind(segs(:,2))>=freq_range(1))
            outliers_h(segs(i,1):segs(i,2)) = false;
            bound(2,1) = max(bound(2,1),ind(segs(i,1)));
            bound(2,2) = min(bound(2,2),ind(segs(i,2)));
        end
        segs = sortrows([ind2seg(outliers_l);ind2seg(outliers_h)],1);
        i_l = 1;	i_r = numel(ind);
        c = c0;
        for i = 1:2
            if any(ind(segs(:,2))<bound(i,1))
                c(1) = i_l == segs(1,1);
                c(3) = ~c(1);
            end
            if any(ind(segs(:,1))>bound(i,2))
                c(2) = i_r == segs(end,2);
                c(4) = ~c(2);
            end
            if c(1) && c(2)
                if i==1
                    c(1) = X(freq_range(1))-x(segs(1,2))>x(segs(end,1))-X(freq_range(2));
                else 
                    c(1) = max(abs(err(segs(1,1):segs(1,2))))>=max(abs(err(segs(end,1):segs(end,2))));
                end
            end
            if i==1	&& any(c(1:4)),	break;	end
        end
        switch find(c,1)
            case 1
                i_l = segs(1,2)+1;
            case 2
                i_r = segs(end,1)-1;
            case 3
                i_l = floor((i_l+segs(1,1))/2)+1;
            case 4
                i_r = ceil((i_r+segs(end,2))/2)-1;
            otherwise
                flag = false;
        end
        ind = ind(i_l:i_r);
        if ~flag,   FIT0;   end
        
        if plt
            for i = 1:4
                if isgraphics(h(i)),	delete(h(i));	end
            end
            it = it+1;
            PLOT;
            title(['Iteration # ',num2str(it)]);
            pause(dt);
%             pause();
        end
    end
else
    x = X;	y = Y;
    FIT0;
    if plt,	PLOT;	end
end

% results
outliers = ind(ind2seg(outliers));
a = exp(P(2));
b = P(1);
fit_range = ind([1,end]);

    function FIT()
        P0 = polyfit(x,y,1);
        outliers0 = y-(P0(1)*x+P0(2))>tDB;
        segs0 = ind2seg(outliers0);
        for j = find(ind(segs0(:,2))<freq_range(1) | ind(segs0(:,1))>freq_range(2))
            outliers0(segs0(j,1):segs0(j,2)) = false;
        end
        P0 = polyfit(x(~outliers0),y(~outliers0),1);
    end

    function FIT0()
        outliers_o = false(size(x));
        for k = 1:length(outliers_o)
            P = polyfit(x(~outliers_o),y(~outliers_o),1);
            outliers = y-(P(1)*x+P(2))>tDB;
            segs0 = ind2seg(outliers);
            for j = find(ind(segs0(:,2))<freq_range(1) | ind(segs0(:,1))>freq_range(2))
                outliers(segs0(j,1):segs0(j,2)) = false;
            end
            if all(outliers==outliers_o)
                break;
            else
                outliers_o = outliers;
            end
        end
    end

    function PLOT()
        figure(fig);
        if p.Results.autofit
            h(1) = plot(x,y,'c.');
            h(2) = plot(x,P0(1)*x+P0(2),'m');
            all_outlier = outliers_l|outliers_h;
            if any(all_outlier)
                h(3) = plot(x(all_outlier),y(all_outlier),'m.');
            end
        end
        if any(outliers)
            plot(x,P(1)*x+P(2),'g');
            plot(x(outliers),y(outliers),'g.');
        end
        if any(outliers0)
            h(4) = plot(x(outliers0),y(outliers0),'r.');
        end
    end

end
