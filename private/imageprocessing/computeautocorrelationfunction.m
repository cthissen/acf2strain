function [Sample] = computeautocorrelationfunction(Sample)
% this function computes the autocorrelation function for a given 2D or 3D
% matrix. The autocorrelation function is shifted such that the center
% value corresponds to lag r = 0.

% calculate estimated autocorrelation function (5-7)
N = numel(Sample.Ts);
scriptT = fftn(Sample.Ts)/N;
scriptR = scriptT.*conj(scriptT);     clear scriptT
Sample.rho = ifftn(scriptR)*N;        clear scriptR
Sample.rho = fftshift(Sample.rho); 

validateattributes(Sample.rho,{'numeric'},{'real'});

end
