%NANSEM Standard error, ignoring NaNs
%   [SE] = NANSEM(X) returns the standard error of X, treating NaNs as
%   missing values. NANSEM operates along the first non-singleton
%   dimension of X. The standard error is computed  by taking the standard
%   deviation with an N-1 normalisation and dividing this by the square
%   root of N, where N is the number of non-NaN elements.
%
%   [M, SE] = NANSEM(X, DIM) takes the standard error along dimension DIM.
%
%   See also NANMEAN, NANSTD.

function [se] = nansem(x, dim)

% Default inputs
if nargin<2
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Take the standard deviation, ignoring nans, and divide by the square root
% of the number of non-nan elements
se = nanstd(x, 0, dim) ./ sqrt(sum(~isnan(x), dim));

end
