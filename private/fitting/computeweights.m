function [w,n]    = computeweights(Sample)
% compute normalized weights for a given weighting scheme

% no weighting
n = numel(Sample.r);
w = ones(size(Sample.r));
       
% normalize weights
w = n * w/sum(w(:));     

%... sanity check
if abs(sum(w)-n) > 1e-6
    warning('sum(w) ~= n ');
    keyboard
end
end