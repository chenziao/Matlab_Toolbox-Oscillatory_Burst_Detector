function y = ncxpdf_cond_gam( x, k, mu, sigma, a, b, prec, Grid, normalize )
% Conditional PDF of noncentral chi distribution given the range of
% non-centrality parameter which follows the Gamma distribution.
% x - values of the chi variable
% k - the shape parameter of the Gamma distribution
% mu - the scale parameter of the Gamma distribution
% sigma - Rayleigh distribution parameter, scaling for the normal variables
% a,b - lower/higher bound of the range of non-centrality parameter [0,Inf]
% prec - precision, the error tolerance for numerical integral
% Grid - preselected grid for interpolating. set to false to disable
% normalize - flag for whether normalize by the probability within bound [a,b]
% y - the conditional probability density values at given x's
if nargin<2 || isempty(k) || k<=0,	k = 2;	end
if nargin<3 || isempty(mu) || mu<=0,	mu = 1;	end
if nargin<4 || isempty(sigma) || sigma<=0,	sigma = 1;	end
if nargin<5 || isempty(a) || a<0,	a = 0;	end
if nargin<6 || isempty(b),	b = Inf;	end
if nargin<7 || isempty(prec),	prec = [];	end
if nargin<8 || isempty(Grid),	Grid = true;	end
if nargin<9 || isempty(normalize),  normalize = true;   end
if a>b,	[b,a] = deal(a,b);	end

% Initialize and find valid points
y = zeros(size(x));
idx = x>0 & x<Inf;

% Special case
if a==b
    X = x(idx)/sigma;
    y(idx) = 2*X.*ncx2pdf(X.^2,2,a^2);	% Noncentral chi distribution
    return;
end

% Normalize parameters
mu = mu/sigma;
a = a/sigma;
b = b/sigma;

% Choose points to evaluate
if isscalar(Grid)
    Interpolate = Grid && numel(x)>300;	% Use interpolation if number larger than 300
    if Interpolate
        kk = max(k,1);
        C = F([a,b],kk+2);	% proportion of power in range
        Mu = (mu^2+2/kk/(kk+1)/(C(2)-C(1)))^0.5;	% Estimate scale
        X = gamma_pdf_int_grid(kk,Mu,prec,max(x(idx))/sigma);
        Grid = X*sigma;
    else
        X = x(idx)/sigma;
    end
else
    X = Grid(:)/sigma;
    Interpolate = true;
end

% % Generate grids for non-centrality parameter
% Z = gamma_pdf_int_grid(k,mu,prec);
% pad = k<1 && a<Z(1);
% Z = [a,Z(Z>a & Z<b)];
% if b<Inf
%     Z = [Z,b];
% else
%     Z = [Z,Z(end)+mu];
% end

% Integral
Y = joint_pb(X);
if normalize
    Y = Y/(F(b,k)-F(a,k));
end

if Interpolate
    y(idx) = interp1(Grid,Y,x(idx),[],0)/sigma;
else
    y(idx) = Y/sigma;
end

    function y = joint_pb(x)
        y = zeros(size(x));
%         % Use trapz on grid
%         fz = gampdf(Z,k,mu);	% pdf of Z
%         if pad
%             fz(1) = fz(2);
%         end
%         for i = 1:numel(x)
%             xz = 2*x(i)*ncx2pdf(x(i)^2,2,Z.^2).*fz;	% joint probability
%             y(i) = trapz(Z,xz);	% integral over [a,b]
%         end
        % Use integral function
        fxz = @(x,z) 2*x*ncx2pdf(x^2,2,z.^2).*gampdf(z,k,mu);	% joint PDF
        for i = 1:numel(x)
            y(i) = integral(@(z) fxz(x(i),z),a,b,'AbsTol',1e-7,'RelTol',1e-4);
        end
    end

    function y = F(z,k)
        y = gammainc(z/mu,k);   % cdf of gamma
    end

end
