function Grid = gamma_pdf_int_grid( k, mu, prec, endpoint )
% Generate grid for numerical intergral of Gamma PDF.
% k - the shape parameter of the Gamma distribution
% mu - the scale parameter of the Gamma distribution
% prec - the error tolerance measured in normalized probability density
% endpoint - the desired maximum value in the grid
if nargin<2 || isempty(mu),	mu = 1;	end
if nargin<3 || isempty(prec),	prec = 0.005;	end

cut_height = 3; % >=2 multiple of tolerance error at cut off point
t = k-1;
gam_k = gamma(k);
prec = min(prec,1/2/cut_height);	% constrain range
decay_rate = (cut_height-1)/cut_height;

% Part 1
[wid1,cut_pt1,tol] = min_binwidth;	% 1st cut off point
if k>=1
    n1 = ceil(cut_pt1/wid1);
    wid1 = cut_pt1/n1;
    grid1 = linspace(0,cut_pt1,n1+1);	% linearly spaced
else
    start_point = cut_pt1*prec;
    grow_rate = gampdf(cut_pt1*decay_rate,k,1);
    n1 = ceil(2*log(gampdf(start_point,k,1))/log(grow_rate))+1;
    grid1 = exp(linspace(log(start_point),log(cut_pt1),n1));
    wid1 = grid1(end)-grid1(end-1);
end
% Part 3
cut_pt2 = inv_pdf(cut_height*tol);	% 2nd cut off point
cut_pt3 = inv_pdf(decay_rate*tol*prec^3);	% end points
n3 = ceil(log(prec^3/cut_height)/log(decay_rate));
wid3 = (cut_pt3-cut_pt2)/n3;
% Part 2
L = cut_pt2-cut_pt1;
if L>wid1 && L>wid3
    r = (L-wid1)/(L-wid3);  % ratio
    n2 = ceil(log(wid3/wid1)/log(r))+1;
    grid2 = cut_pt1+cumsum(cumprod([wid1,r*ones(1,n2-1)])); % geometric spaced
    cut_pt2 = grid2(end);
else
    n2 = 0; grid2 = [];
    cut_pt2 = cut_pt1;
end
% Part 3
grid3 = cut_pt2+(1:n3)*wid3;	% linearly spaced
cut_pt3 = grid3(end);

Grid = [grid1,grid2,grid3]*mu;
if nargin>=4
    idx = find(Grid<endpoint,1,'last');
    Grid(idx+1:end) = [];
    Grid = [Grid,endpoint];
else
    Grid(~isfinite(Grid)) = [];
end

% disp([1;1/k]*[cut_pt1,cut_pt2,cut_pt3]);
% disp([n1,n2,n3]);

    % Pick minimum bin width depending on the maximum slope of gamma pdf.
    function [binwidth,min_slp_x,tol] = min_binwidth()
        if k>=1
            % Find maximum decaying slope
            min_slp_x = t+t^0.5;
            max_slp = gampdf(min_slp_x,k,1)/(1+t^0.5);	% slope = -max_slp
            % Want: bin width * min(slope) <= maximum pdf * precision
            max_pdf = (t*exp(-1))^t/gam_k;	% maximum pdf depending on k
            tol = max_pdf*prec;     % error tolerance in probability density
            binwidth = tol/max_slp;
        else
            min_slp_x = inv_pdf(1);
            max_slp = 1-t/min_slp_x;
            tol = prec;
            binwidth = tol/max_slp;
        end
    end

    % Inverse function of gampdf with mu=1 in the tail
    function x = inv_pdf(y)
        if k == 1
            x = -log(y);
        else
            branch = -1*(k>1);
            try
                x = -t*real(lambertw(branch,-(y*gam_k).^(1/t)/t));
                flag = ~isfinite(x);
            catch
                flag = true;
            end
            if flag
                x = fminsearch(@(x) err_inv_pdf(x,y),k,optimset('TolFun',0.01*y));
            end
        end
    end

    function err = err_inv_pdf(x,y)
        if x>=t
            err = abs(gampdf(x,k,1)-y);
        else
            err = 1;
        end
    end
end
