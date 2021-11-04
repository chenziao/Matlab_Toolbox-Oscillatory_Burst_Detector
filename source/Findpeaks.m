function [ peaks, idx ] = Findpeaks( x, threshold )
trend = [sign(diff(x(:)));1];
f_idx = find(trend==0);
for i = length(f_idx):-1:1,
    if trend(f_idx(i)+1)>=0, trend(f_idx(i)) = 1;	else trend(f_idx(i)) = -1; end
end
if nargout > 1
    idx = find([0;diff(trend)]==-2);
    peaks = x(idx);
    if nargin > 1
        idx = idx(peaks>=threshold);
        peaks = x(idx);
    end
else
    peaks = x([0;diff(trend)]==-2);
    if nargin >1
        peaks = peaks(peaks>=threshold);
    end
end
end
