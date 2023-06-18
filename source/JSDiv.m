function div = JSDiv(P,Q,alpha)
% div = JSD(P,Q) Jensen-Shannon divergence of two probability distributions
% P and Q are automatically normalised to have the sum of one on rows
% P, Q are vectors with same length
% alpha is a value between 0 and 1 added to the count for each event,
% so called add-alpha smoothing, also called Laplace smoothing. default:0
% If 0 count is found, alpha is automatically 1 if not specified.
if nargin<3 || ~isscalar(alpha) || alpha<=0 || alpha>1
    alpha = 1-all(P|Q);
end
P = P+alpha;
Q = Q+alpha;
P = P/sum(P);
Q = Q/sum(Q);
M = (P+Q)/2;
div = (KLDiv(P,M,'JS')+KLDiv(Q,M,'JS'))/2;
end