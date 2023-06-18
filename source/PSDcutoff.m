function PSD = PSDcutoff( PSD,range,level )
% Calculate cutoff line of PSD (with linear frequency points including 0 Hz)
% PSD - original power spetral density
% range - (n-by-2) index of segment range for cutoff
% level - cutoff level, same size with PSD, inferred linearly from boundary points if not specified.
% Return PSD after resetting power to cutoff level over the range
if nargin<3
    m = numel(PSD);
    intvl = range+[-1,1];
    intvl = [max(intvl(:,1),2),min(intvl(:,2),m)];
    range = [max(range(:,1),2),min(range(:,2),m)];
    level = reshape(PSD(intvl),size(intvl));
    for i = 1:size(range,1)
        idx = range(i,1):range(i,2);
        if intvl(i,1)==intvl(i,2)
            PSD(idx) = level(i,1);
        else
            PSD(idx) = exp(interp1(log(intvl(i,:)-1),log(level(i,:)),log(idx-1)));
        end
    end
else
    for i = 1:size(range,1)
        idx = range(i,1):range(i,2);
        PSD(idx) = level(idx);
    end
end
end