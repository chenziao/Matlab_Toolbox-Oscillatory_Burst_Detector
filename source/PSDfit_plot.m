function ax = PSDfit_plot(results,flag,ax,zoom_in)
if nargin<2 || isempty(flag),	flag = 1;	end
if nargin<3,	ax = [];	end
if nargin<4 || isempty(zoom_in),	zoom_in = false;	end

f = results.f;
PSD = results.PSD;
PSD_smoo = results.PSD_smoo;
fit_i = results.fit_i;
PSD_fit = results.fit_a*f.^results.fit_b;
outliers = results.outliers;
nDB = results.nDB;
i_fit_range = results.i_fit_range;
freq_range = results.freq_range;

if isempty(ax),	figure;	ax = gca;	end
hold(ax,'on');	h = zeros(1,6);
h(1) = plot(ax,f,pow2db(PSD),'b');
h(2) = plot(ax,f,pow2db(PSD_smoo),'m');
if flag
    h(4) = plot(ax,f(fit_i),pow2db(PSD_fit(fit_i)),'c--','LineWidth',2);
    for i = 1:size(outliers,1)
        o_i = outliers(i,1)-1:outliers(i,2)+1;
%         cutoff = PSDcutoff(f(o_i),PSD_smoo(o_i));
        cutoff = PSDcutoff(f(o_i),[],[],db2pow(nDB)*PSD_fit(outliers(i,:))');
        h(5) = plot(ax,f(outliers(i,1):outliers(i,2)),pow2db(cutoff(2:end-1)),'g--','LineWidth',2);
    end
end
if flag == 1
    h(3) = plot(ax,f(i_fit_range(1):fit_i(1)-1),pow2db(PSD_fit(i_fit_range(1):fit_i(1)-1)),'r--','LineWidth',2);
    h(3) = plot(ax,f(fit_i(end)+1:i_fit_range(2)),pow2db(PSD_fit(fit_i(end)+1:i_fit_range(2))),'r--','LineWidth',2);
else
    h(3) = plot(ax,f(i_fit_range(1):i_fit_range(2)),pow2db(PSD_fit(i_fit_range(1):i_fit_range(2))),'r--','LineWidth',2);
end
[~,f_i] = max(PSD_smoo);    xl = [min(f(f_i),2),f(end-1)];
if zoom_in
    xl(1) = max(f(i_fit_range(1)),xl(1));
    xl(2) = min(f(i_fit_range(2)),xl(2));
end
xlim(ax,xl);    ylim(ax,'auto');	yl = get(ax,'YLim');
h(6) = plot(ax,freq_range(1)*[1,1],yl,'y:');	plot(ax,freq_range(2)*[1,1],yl,'y:');
set(ax,'XScale','log');
leg = {'raw PSD','smoothed PSD','fit curve','fit background component',...
    'signal frequency range','frequency bound for signal of interest'};
if ~flag
    h = h([1:3,6]);	leg = leg([1:3,6]);
end
set(legend(ax,h,leg,'Location','southwest'),'Box','off','Color','None');
xlabel(ax,'Freguency (Hz)');	ylabel(ax,'PSD (dB/Hz)');
hold(ax,'off');

end