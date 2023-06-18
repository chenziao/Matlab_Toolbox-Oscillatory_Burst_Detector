function results = getBurstAttr(results,LFP,fs,freq_range,Nfft)
%% Find frequency
tic;
duridx = results.duridx;
duration = double(duridx(:,2)-duridx(:,1)+1);
nmerged = results.nmerged;
burst_freq = zeros(nmerged,1);
cyc_num = zeros(nmerged,1);

freq = (0:Nfft/2)/Nfft*fs;
fid = max(find(freq>=freq_range(1),1,'first')-1,1):min(find(freq<=freq_range(2),1,'last')+1,Nfft/2);
for i = 1:nmerged
    if duridx(i,2)-duridx(i,1)>Nfft-3
        idx = max(results.tpk(i)-Nfft/2,duridx(i,1));
        idx = min(idx,duridx(i,2)-Nfft+3);
        x = LFP(idx:idx+Nfft-3);
    else
        x = LFP(duridx(i,1):duridx(i,2));
    end
    x = x-mean(x);
    x_fft = fft([0;x;0].*tukeywin(numel(x)+2,0.25),Nfft);
    [M,I] = Findpeaks(abs(x_fft(fid)));
    if ~isempty(I)
        [~,Im] = max(M);
        burst_freq(i) = freq(fid(I(Im)));
        cyc_num(i) = duration(i)*burst_freq(i)/fs;
    end
end
ivalid = burst_freq>0;
nvalid = sum(ivalid);
runtime2 = toc;

%% Results
results = rmfield(results,'tpk');
results.nvalid = nvalid;
results.bursts_seg_id = int32(zeros(nvalid,1));
results.duridx = duridx(ivalid,:);
results.AP = results.AP(ivalid);
results.CN = cyc_num(ivalid);
results.BF = burst_freq(ivalid);
results.df = fs/Nfft;
results.runtime2 = runtime2;

end