function param = fitdist_trunc( data,pdffun,cdffun,min_max,start )
% Fitting the truncated customized distribution to the data
% data - input data
% pdffun - pdf function of desired distribution type
% cdffun - pdf function of desired distribution type
% min_max - range verctor defining the border where the data is truncated
% (Use -+Inf for the side(s) not truncated, NaN for that defined by the data limit.)
% start - vector of initial parameter values for maximum likelihood fitting

if isempty(min_max)
    min_max = [min(data),max(data)];  
else
    if numel(min_max)~=2
        error('The range should contain both min and max');
    end
    if isnan(min_max(1))
        min_max(1) = min(data);  
    end
    if isnan(min_max(2))
        min_max(2) = max(data);  
    end
end
ind = data>=min_max(1) & data<=min_max(2);
if sum(ind)==0
    error('The fitting range [%g,%g] is empty',min_max(1),min_max(2));
else
    data = data(ind);
end

% find the maximum likelihood estimate
param = mle(data,'pdf',@pdf_trunc,'start',start);

    function y = pdf_trunc(x,varargin)
        y = pdffun(x,varargin{:})./(cdffun(min_max(2),varargin{:})- ...
            cdffun(min_max(1),varargin{:}));
    end

end
