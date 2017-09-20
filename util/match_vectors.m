function idx = match_vectors(a,b)

% idx = MATCH_VECTORS(a,b) - For two vectors a and b in which b is a 
%   permutation of a, match_vectors returns the indices idx such that
%   a = b(idx).  

if ~isequal(class(a), class(b))
    error('Input vectors are not of the same class');
end

[~, aix] = sort(a);
[~, bix] = sort(b);

aix_inv(aix) = 1:numel(a);
idx = bix(aix_inv);

end