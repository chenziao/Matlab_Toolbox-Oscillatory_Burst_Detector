function results = getDurations(LFP_amp,threshold,StopProportion)
tic;
%% detect peaks
T = numel(LFP_amp);
[peaks,tpks,trend] = Findpeaks(LFP_amp);
tpks = int32(tpks);
npeak = numel(peaks);
ipk = find(peaks>threshold);
peak = peaks(ipk);
tpk = tpks(ipk);
ndetect = numel(tpk);

%% Find duration algorithm
% Search time points along time direction. Store the value of a peak that
% the search has passed into a sorted list. Compare value at each point
% with the certain cutoff value of the maximum peak in the list. Remove
% from the list when the cutoff edge of a peak is found. Revert the search
% time direction for finding the edge on the opposite side.
ipk = [0;ipk;numel(tpks)+1];
tpksx = [0;tpks;T+1];
trend = [-1;trend;1];
sidestops = StopProportion*peak;
nextid = zeros(ndetect+1,1);
head = ndetect+1;
duridx = zeros(ndetect,2);

% forward path
for i = 1:ndetect
    pre = head; ptr = nextid(head);
    while ptr>0 && sidestops(i)<sidestops(ptr)
        pre = ptr;
        ptr = nextid(pre);
    end
    nextid(i) = nextid(pre);
    nextid(pre) = i;
    j = ipk(i+1);
    while j<ipk(i+2)
        t = tpksx(j+1)+1;
        while t<tpksx(j+2)
            while LFP_amp(t)<=sidestops(nextid(head))
                duridx(nextid(head),2) = t-1;
                nextid(head) = nextid(nextid(head));
                if nextid(head)==0
                    j = ipk(i+2)-1;	t = tpksx(j+2)-1;	break;
                end
            end
            t = t+1;
            if trend(t)>0,	break;  end
        end
        j = j+1;
    end
end
while nextid(head)>0
    duridx(nextid(head),2) = T;
    nextid(head) = nextid(nextid(head));
end

% backward path
for i = ndetect:-1:1
    pre = head; ptr = nextid(head);
    while ptr>0 && sidestops(i)<sidestops(ptr)
        pre = ptr;
        ptr = nextid(pre);
    end
    nextid(i) = nextid(pre);
    nextid(pre) = i;
    j = ipk(i+1)+1;
    while j>ipk(i)+1
        t = tpksx(j)-1;
        while t>tpksx(j-1)
            while LFP_amp(t)<=sidestops(nextid(head))
                duridx(nextid(head),1) = t+1;
                nextid(head) = nextid(nextid(head));
                if nextid(head)==0
                    j = ipk(i)+2;	t = tpksx(j-1)+1;	break;
                end
            end
            if trend(t)<0,	break;  end
            t = t-1;
        end
        j = j-1;
    end
end
while nextid(head)>0
    duridx(nextid(head),1) = 1;
    nextid(head) = nextid(nextid(head));
end

%% Merge overlapping duration algorithm
% For any two distinct peaks, there are only two possible relations between
% their duration intervals, disjoint (non-ovelapping) or inclusion (overlapping).
% The main idea is to remove the smaller one if two peaks have overlapping durations.
ivalid = true(ndetect,1);  % valid peaks
i = 1;  j = 1;
while j<ndetect
    j = j+1;	% next j with peak time to the right
    if duridx(i,2)<duridx(j,1)
        i = j;	% duration i disjoint with j, update i to current j for next comparison
    elseif peak(i)<peak(j)
        ivalid(i) = false;	% duration i includes j, keep j
        i = j;	% update i to current j for next comparison
    else
        ivalid(j) = false;	% duration j includes i, keep i
    end
end
duridx = int32(duridx(ivalid,:));
runtime1 = toc;

%% Results
results = struct('peaks',peaks,'tpks',tpks,'pkpha',[],'pks_seg_id',int32(zeros(npeak,1)), ...
    'npeak',npeak,'ndetect',ndetect,'nmerged',size(duridx,1),'tpk',tpk(ivalid),'nvalid',[], ...
    'bursts_seg_id',[],'duridx',duridx,'AP',peak(ivalid),'CN',[],'BF',[],'df',[]);
results.runtime1 = runtime1;

end