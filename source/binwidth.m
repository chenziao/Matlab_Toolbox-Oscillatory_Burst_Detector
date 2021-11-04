function d = binwidth( x )
% Calculate the average bin width given bin centers or edges in vector x.
if ~isvector(x) || numel(x)<2
    error('Input should be a vector.');
end
d = (x(end)-x(1))/(length(x)-1);
end

