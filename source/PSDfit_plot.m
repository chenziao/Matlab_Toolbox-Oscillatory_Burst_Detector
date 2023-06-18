function ax = PSDfit_plot(results,flag,ax,zoom_in)
if nargin<2 || isempty(flag),	flag = 1;	end
if nargin<3,	ax = [];	end
if nargin<4 || isempty(zoom_in),	zoom_in = false;	end

f = results.f;
PSD = results.PSD;
PSD_smoo = results.PSD_smoo;
bg_range_i = results.bg_range_i;
PSD_fit = getBackgroundFit(results);
outliers = results.outliers;
tDB = results.tDB;
i_fit_range = results.i_fit_range;
freq_band = results.freq_band;

if isempty(ax),	figure;	ax = gca;	end
hold(ax,'on');	h = zeros(1,6);
h(1) = plot(ax,f,pow2db(PSD),'b');
h(2) = plot(ax,f,pow2db(PSD_smoo),'m');
if flag
    h(4) = plot(ax,f(bg_range_i),pow2db(PSD_fit(bg_range_i)),'c--','LineWidth',2);
    for i = 1:size(outliers,1)
        o_i = outliers(i,1):outliers(i,2);
        cutoff = pow2db(PSD_fit(o_i))+tDB;
        h(5) = plot(ax,f(o_i),cutoff,'g--','LineWidth',2);
    end
end
if flag == 1
    h(3) = plot(ax,f(i_fit_range(1):bg_range_i(1)), ...
        pow2db(PSD_fit(i_fit_range(1):bg_range_i(1))),'r--','LineWidth',2);
    h(3) = plot(ax,f(bg_range_i(end):i_fit_range(2)), ...
        pow2db(PSD_fit(bg_range_i(end):i_fit_range(2))),'r--','LineWidth',2);
else
    h(3) = plot(ax,f(i_fit_range(1):i_fit_range(2)), ...
        pow2db(PSD_fit(i_fit_range(1):i_fit_range(2))),'r--','LineWidth',2);
end
xl = [f(2),f(end-1)];
if zoom_in
    xl(1) = max(f(i_fit_range(1)),xl(1));
    xl(2) = min(f(i_fit_range(2)),xl(2));
end
xlim(ax,xl);    ylim(ax,'auto');	yl = get(ax,'YLim');
h(6) = plot(ax,freq_band(1)*[1,1],yl,'--','color','#EDB120');
plot(ax,freq_band(2)*[1,1],yl,'--','color','#EDB120');
set(ax,'XScale','log');
leg = {'Raw PSD','Smoothed PSD','Fit curve','Fit background component',...
    'Signal frequency range','Signal frequency band'};
if ~flag
    h = h([1:3,6]);	leg = leg([1:3,6]);
end
set(legend(ax,h,leg,'Location','southwest'),'Box','off','Color','None');
xlabel(ax,'Frequency (Hz)');	ylabel(ax,'PSD (dB/Hz)');
hold(ax,'off');

end