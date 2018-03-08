function GNobj = GNgenerate(t, N, varargin)

if nargin > 2 && ~isempty(varargin{1})
    GNmean = varargin{1};
else
    GNmean = 0;
end
if nargin > 3 && ~isempty(varargin{2})
    GNvar = varargin{2};
else
    GNvar = 1;
end

if length(GNmean)~= length(GNvar)
    error('Sizes of mean and var are not equal');
else
    GNnoise = randn(t, N); 
    GNnoise = GNnoise - repmat(mean(GNnoise), t, 1);
    GNnoise = GNnoise.*repmat(1./(std(GNnoise)), t,1);
    if length(GNmean) <=1
        GNnoise = repmat(GNmean, t, N)...
            + repmat(sqrt(GNvar), t, N).*GNnoise;
    else
        GNnoise = repmat(GNmean, t,1)...
            + repmat(sqrt(GNvar), t, 1).*GNnoise;
    end
    GNobj.GNnoise = GNnoise;
    GNobj.UPlimit = max(GNnoise);
    GNobj.LOWERlimit = min(GNnoise);
end

end

