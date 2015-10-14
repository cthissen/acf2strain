function [zeta0,c]          = computezeta0(r0,Sample)
% compute best-fit for isotropic autocorrelation function

% splinefit is a program available from the Matlab user area.
% This program is based on a user-specified set of knots (x positions), 
% and finds the spline function using a least-squares criterion.

breaks = linspace(min(r0),max(r0),Sample.knots);
[c]    = splinefit(r0,Sample.zeta,breaks);
zeta0  = ppval(c,r0);
        
end