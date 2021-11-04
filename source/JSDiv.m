function dist = JSDiv(P,Q,alpha)
% dist = JSD(P,Q) Jensen-Shannon divergence of two probability distributions
% P and Q are automatically normalised to have the sum of one on rows
% P, Q are vectors with same length
% alpha is a value between 0 and 1 added to the count for each event,
% so called add-alpha smoothing, also called Laplace smoothing. default:0
% If 0 count is found, alpha is automatically 1 if not specified.
if ~isvector(P)||~isvector(Q)
    error('P and Q should be vectors.');
end
if length(P)~=length(Q)
    error('The length of P and Q should be the same.');
end
if any(P<0)||any(Q<0)||~all(isfinite(P))||~all(isfinite(Q))
    error('The inputs should be finite non-negative values.');
end
if nargin<3 || ~isscalar(alpha) || alpha<=0 || alpha>1
    alpha = 1-all(P|Q);
end

Q = Q+alpha;
P = P+alpha;
Q = Q/sum(Q);
P = P/sum(P);
M = (P+Q)/2;
dist = (KLDiv(P,M,'JS')+KLDiv(Q,M,'JS'))/2;
end