function div = KLDiv(P,Q,alpha)
% div = D_KL(P||Q) Kullback-Leibler divergence of two discrete probability distributions.
% P and Q are automatically normalised to have the sum of one on rows
% P, Q are vectors of counts with same length
% alpha is a value between 0 and 1 added to the count for each event,
% so called add-alpha smoothing, also called Laplace smoothing. default:0
% If 0 count is found, alpha is automatically 1 if not specified.
% alpha='JS' is used only for efficiency when calculating Jensen-Shannon divergence.
if nargin<3 || ~strcmp(alpha,'JS')
    if nargin<3 || ~isscalar(alpha) || alpha<=0 || alpha>1
        alpha = 1-all(P&Q);
    end
    P = P+alpha;
    Q = Q+alpha;
    P = P/sum(P);
    Q = Q/sum(Q);
end
p = P./Q;
div = sum(P.*log(p));
end