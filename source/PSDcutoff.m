function PSD = PSDcutoff( f,PSD,range,level )
% Calculate cutoff line of PSD which is a straight line in log-log scale
% f - frequency points
% PSD - power spetral density
% range - (n-by-2) index of range in f to cut (default: full range [2,end-1])
% level - if specified, corresponds to cutoff level at boundary points of range
% Return PSD after removing power above cutoff line
if nargin<2 || isempty(PSD)
    PSD = zeros(size(f));
end
if nargin<3 || isempty(range)
    range = [2,numel(f)-1];
end
if nargin<4 || size(level,2)~=2
    intvl = bsxfun(@plus,range,[-1,1]);
    level = reshape(PSD(intvl),size(intvl));
else
    intvl = range;
end
lgf = reshape(log(f),1,[]);
for i = 1:size(range,1)
    idx = range(i,1):range(i,2);
    if intvl(i,1)==intvl(i,2)
        PSD(idx) = level(i,1);
    else
        PSD(idx) = exp(interp1(lgf(intvl(i,:)),log(level(i,:)),lgf(idx)));
    end
end
end