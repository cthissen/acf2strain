function [Sample] = reducedata(Sample)
%% Reduce data
% reduce data by removing symmetric part of ACF and by eliminating the r=0
% lag

Sample.rx = Sample.rx(:);
Sample.ry = Sample.ry(:);
Sample.rz = Sample.rz(:);

switch Sample.d
    case 2
        uniqueIdx = find(Sample.ry < 0);        
    case 3
        uniqueIdx = find(Sample.rz < 0);
end
Sample.zeta(uniqueIdx) = [];
Sample.r(uniqueIdx) = [];
Sample.rx(uniqueIdx) = [];
Sample.ry(uniqueIdx) = [];
Sample.rz(uniqueIdx) = [];


% reduce data by removing large lags for faster fitting
largeLagIdx = find(Sample.r(:) > Sample.maxRadialLag);
Sample.zeta(largeLagIdx) = [];
Sample.r(largeLagIdx) = [];
Sample.rx(largeLagIdx) = [];
Sample.ry(largeLagIdx) = [];
Sample.rz(largeLagIdx) = [];

% also remove r = 0 (rho ~= 1, z = inf)
infIdx = find(Sample.r == 0);
Sample.zeta(infIdx) = [];
Sample.r(infIdx) = [];
Sample.rx(infIdx) = [];
Sample.ry(infIdx) = [];
Sample.rz(infIdx) = [];

end