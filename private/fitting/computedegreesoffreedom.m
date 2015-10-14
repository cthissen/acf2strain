function [nu] = computedegreesoffreedom(Sample,n)
% compute degrees of freedom given weighting scheme and dimensionality

switch Sample.d
    case 2
        mE = 2; % using isochoric constraint
    case 3
        mE = 5; % using isochoric constraint
end



mC = 4*(Sample.knots);

nu = n-mE-mC;

end
