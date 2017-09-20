function rnk = rankorder(m, varargin)
% RANKORDER Compute rank order of elements in a numeric matrix.
%   RNK = RANKORDER(M)
%   Ranks are computed by sorting M in ascending order along each column.
%   RANKORDER(..., 'param1', 'value1',...) specify optional parameters.
%   Valid parameters are:
%
%   'as_fraction'   Logical, returns ranks as a fraction of total rows (or
%   columns
%   'as_percentile'   Logical, returns ranks as a percentile of total rows (or
%   columns
%   'dim'           Integer, sorts along dimension specified can be 1 for
%                   rows or 2 for columns (is 1 by default).
%   'direc'         String, sort in ascending or descending order. Can be
%                   'ascend' or 'descend' (is 'ascend' by default).
%   'fixties'       Logical, adjusts for ties (is true by default).
%   'zeroindex'     Logical, if true the returned ranks are zero indexed
%                   (is false by default).

% $Author: Rajiv Narayan [narayan@broadinstitute.org]
% $Date: Jul.01.2010 12:01:46 EDT

% Changes: wrapper for tiedrank now

% nin=nargin;
% error(nargchk(1,3,nin,'struct'));

% breakties = false;

pnames = {'dim', 'direc', 'zeroindex',...
          'fixties', 'as_fraction', 'as_percentile'};
dflts = {1, 'ascend', false,...
        true, false, false};
arg = parse_args(pnames, dflts, varargin{:});

if arg.fixties
    rank_alg = @tiedrank;
else
    rank_alg = @myrankorder;
end

if arg.as_fraction || arg.as_percentile
    arg.zeroindex = true;
end

if ndims(m)<=2
    if isequal(arg.dim, 2)
        m = m';
    end
    switch lower(arg.direc)
        case 'ascend'
            rnk = rank_alg(m);
        case 'descend'
            rnk = rank_alg(-m);
        otherwise
            error ('Invalid direc specified: %s', arg.direc);
    end
    
    if arg.zeroindex
        rnk = rnk-1;
    end
    if arg.as_fraction || arg.as_percentile
       nrows = size(rnk, 1)-1;
       scale = (1 + 99*arg.as_percentile);
       if nrows > 0
           rnk = scale*rnk / nrows;
       else
           rnk = scale;
       end       
    end
    if isequal(arg.dim, 2)
        rnk = rnk';
    end
else
    error('Input should have <=2 dimensions')
end
end

% Compute ranks for matrix m along the first dimension
% Does not adjust for ties
function m = myrankorder(m)
[nr, nc]= size(m);
if nr==1 && nc>1 
    m=m(:);
end
ord = (1:nr)';
[~, sidx] = sort(m, 1);
sidx = bsxfun(@plus, sidx, (nr*(0:nc-1)));
% m = zeros(nr, nc);
% m(sidx) = ord(:, ones(nc,1));
m(sidx) = repmat(ord, 1, nc);
end
