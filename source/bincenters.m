function ctrs = bincenters( x, add )
% Calculate the center of bins given bin edges in evenly spaced vector x.
% Insert points in vector add.
if ~isvector(x) || numel(x)<2
    error('Input should be a vector.');
end
ctrs = x(1:end-1)+binwidth(x)/2;
if nargin>1
    ctrs = sort([add, ctrs]);
end
