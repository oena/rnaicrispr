function vec = mat2vec(mymat, varargin)

% VEC = MAT2VEC(MYMAT, VARARGIN)
%   Takes as input an nxm numerical matrix mymat and returns a vector vec
%   of size (nxm)x1.  
%   
%   Arguments:
%     usecol - boolean, 1 to return column vector; 0 to return row
%         vector.  Default value is 1.
%     rmdiag - boolean, 1 to omit the diagonal, 0 to return the full matrix.  Default 0.
%         Only applies if the matrix is square, otherwise rmdiag is not applied.


pnames = {'usecol', 'rmdiag'};
dflts = {1, 0};
args = parse_args(pnames, dflts, varargin{:});

if and(args.rmdiag, (size(mymat,1) == size(mymat,2)))
  vec = mymat; 
  vec(1:size(vec,1)+1:end) = [];
else
  vec = reshape(mymat, numel(mymat), 1);
end


if ~args.usecol
    vec = vec';
end

end
