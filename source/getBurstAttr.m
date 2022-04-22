function results = getBurstAttr(LFP,LFP_ASamp,fs,threshold,freq_range,StopProportion,Nfft)
T = numel(LFP);

% detect peaks
[peaks,tpks] = Findpeaks(LFP_ASamp);
tpks = int32(tpks);
tpk = tpks(LFP_ASamp(tpks)>threshold);
peak = LFP_ASamp(tpk);
ndetect = numel(tpk);

% Find duration
duridx = zeros(ndetect,2);
for i = 1:ndetect
    sidestop = StopProportion*peak(i);
    t = tpk(i)-1;
    while t>0&&LFP_ASamp(t)>sidestop, t = t-1; end
    duridx(i,1) = t+1;
    t = tpk(i)+1;
    while t<=T&&LFP_ASamp(t)>sidestop, t = t+1; end
    duridx(i,2) = t-1;
end

% Merge overlap duration
d = diff([false;duridx(1:end-1,2)>duridx(2:end,1);false]);
td = [find(d==1),find(d==-1)];
vpkid = true(ndetect,1);
for i = 1:size(td,1)
    id = td(i,1):td(i,2);
    np = numel(id);
    [~,I] = sort(peak(id),'descend');	% sort peaks by peak value
    id = id(I);
    t = tpk(id); % time of sorted peaks
    for j = 1:np
        Id = id(j);
        if vpkid(Id)
            idx = j+1:np;   % index for smaller peaks
            vpkid(id(idx((t(idx)>t(j)&duridx(id(idx),1)<duridx(Id,2))|...
                (t(idx)<t(j)&duridx(id(idx),2)>duridx(Id,1))))) = false;
        end
    end
end
duridx = duridx(vpkid,:);
duration = duridx(:,2)-duridx(:,1)+1;
nmerged = size(duridx,1);

% Find frequency
freq = (0:Nfft/2)/Nfft*fs;
fid = find(freq>=freq_range(1)&freq<=freq_range(2));
burst_freq = zeros(nmerged,1);
cycnum = zeros(nmerged,1);
for i = 1:nmerged
    y = LFP(duridx(i,1):duridx(i,2));
    if numel(y)>Nfft-2
        idx = floor((numel(y)-Nfft+2)/2);
        y = y(idx+1:idx+Nfft-2);
    end
    y = y-mean(y);
    yfft = abs(fft([0;y;0].*tukeywin(numel(y)+2,0.25),Nfft));
    [~,I] = Findpeaks(yfft(fid));
    if ~isempty(I)
        [~,Im] = max(yfft(fid(I)));
        burst_freq(i) = freq(fid(I(Im)));
        cycnum(i) = duration(i)*burst_freq(i)/fs;
    end
end

% Results
npeak = numel(peaks);
peak = peak(vpkid);
ivalid = burst_freq>0;
nvalid = sum(ivalid);
results = struct('peaks',peaks,'tpks',tpks,'npeak',npeak,'ndetect',ndetect,'nmerged',nmerged,...
    'nvalid',nvalid,'AP',peak(ivalid),'CN',cycnum(ivalid),'BF',burst_freq(ivalid),'df',fs/Nfft);

end