function [Sample] = estimategrainsize(Sample)
% This function estimates the grain size for a given istropic ACF. Three
% models are implemented, an exponential decay, r^*, and R_S. See Thissen
% and Brandon (2015) for further details.

%% Notes
%... Sequential fitting for observed ACF:
% 1) Spherical model
%   fit zeta = real( atanh( DRho ...
%     + (1-DRho).*(r<=2*Rs).*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3) ));
%   start and constraints:
%      DRho:  DRhoEst, [-0.999, + 0.999]
%      Rs:    Rmean04, [+1e-6, +Inf]

% 2) Two-Component model: A modified version with
% the spherical model for the self-correlation component, and
% a decaying sinusoid for the neighbor-correlation component.
%   fit zeta = real( atanh( DRho ...
%      + (1-DRho).*((r<=2*Rs).*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3) ...
%      + A*exp(-(r/Rn).^Beta)*sin(pi*(r/Rn)))));
%   start and constraints:
%      A:     0, [-Inf, +Inf]
%      Beta:  1, [-Inf, +Inf]
%      DRho:  DRhoEst, [-0.999, + 0.999]
%      Rn:    5*Rmean04, [+1e-6, +Inf]
%      Rs:    5*Rmean04, [+1e-6, +Inf]


%% Initialize the data and some preliminary estimates
r = Sample.r0;
rho = tanh(Sample.zeta);
fid = [Sample.fid];

%... Make sure that data are sorted with increasing r
r = r(:);
rho = rho(:);
[r, k] = sort(r);
rho = rho(k);

%... Trim datat to remove entries that have r==0 and Rho == 1, 
% which causes problems with the zeta transform
i = (r~=0) & (rho<1);
r = r(i);
rho = rho(i);
zeta = atanh(rho);

dualfprintf(fid,'\n\n======= Grain size estimate =======\n')

%% Fit using approximation from Panozzo Heilbronner (1992)

%... approximate ACF using truncated Taylor Series
k = find(rho<0, 1) - 1;
% Define (rho-1) = G*coeffRmean
G = [r(1:k), r(1:k).^2, r(1:k).^3, r(1:k).^4, r(1:k).^5, r(1:k).^6];
coeffRmean = G\(rho(1:k) - 1);
rPredRmean = linspace(0,r(k))';
G = [rPredRmean, rPredRmean.^2, rPredRmean.^3, rPredRmean.^4, rPredRmean.^5, rPredRmean.^6];
rhoPredRmean = G*coeffRmean + 1;

rStar = interp1(rhoPredRmean, rPredRmean, 1/2);
Rmean04 = 1.44*rStar;
% dualfprintf(fid,'Estimate Rmean using approximation from Panozzo Heilbronner (1992)\n');
% dualfprintf(fid,sprintf('Rmean = %.3f\n\n', Rmean04));



%% Fit using Spherical Model
modelName = 'Spherical Model';
% Model function:
% zeta = real( atanh( DRho + ...
% (1-DRho)*(r<=2*Rs)*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3) ));
ft = fittype( ...
    ['real( atanh( + ', ...
    '(r<2*Rs)*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3) ))'], ...
    'independent', 'r', 'dependent', 'zeta' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%... Weights
opts.Weights = ones(size(r));
%... Set starting point for parameter search.
% Order: DRho, Rs
opts.StartPoint = [Rmean04];
%... Set constraints for ranges for best-fit parameters.
opts.Lower = [1e-6];
opts.Upper = [+Inf];

% Fit model to data
[fitResultSM, gofSM] = fit( r, zeta, ft, opts);
coeff = coeffvalues(fitResultSM);
names = coeffnames(fitResultSM);
dualfprintf(fid,sprintf('BEST-FIT SOLUTION: %s\n', modelName));
% dualfprintf(fid,sprintf('Filename: %s\n', fileName));
for i = 1:length(coeff)
    dualfprintf(fid,sprintf('%8s: %10.3f\n',char(names(i)),coeff(i)));
end
dualfprintf(fid,sprintf('%8s: %10.3f\n\n','R^2',  gofSM.rsquare));

%% Fit using Two-Component Model
modelName = 'Two-Component model';
% Model function:
% zeta = real( atanh( DRho ...
% + (1-DRho).*((r<=2*Rs).*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3) ...
% + fn*exp(-(r/Rn).^Beta).*sin(pi*r/Rn)) ));
ft = fittype( ...
    ['real( atanh( ', ...
    '+ ((r<=2*Rs).*(1 - (3/4)*r/Rs + (1/16)*(r/Rs).^3)', ...
    '+ fn*exp(-(r/Rn).^Beta).*sin(pi*r/Rn)) ))'], ...
    'independent', 'r', 'dependent', 'zeta' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%... Weights
opts.Weights = ones(size(r));
%... Set starting point for parameter search.
% Order: Beta, DRho, Rn, Rs, fn
opts.StartPoint = [1 5*Rmean04 5*Rmean04 0];
%... Set constraints for ranges for best-fit parameters
opts.Lower = [-Inf +1e-6 +1e-6 -Inf];
opts.Upper = [+Inf +Inf  +50  +Inf];

% Fit model to data.
[fitResultTCM, gofTCM] = fit( r, zeta, ft, opts);
coeff = coeffvalues(fitResultTCM);
names = coeffnames(fitResultTCM);
dualfprintf(fid,sprintf('BEST-FIT SOLUTION: %s\n', modelName));
for i = 1:length(coeff)
    dualfprintf(fid,sprintf('%8s: %10.3f\n',char(names(i)),coeff(i)));
end
dualfprintf(fid,sprintf('%8s: %10.3f\n\n','R^2', gofTCM.rsquare));

%... Test to find out whether the two-component model provides
% a significant improvement in fit, using the F test from
% p. 207 in Bevington and Robinson, 2003.
df1 = 1;
df2 = length(zeta) - 5;
F = (gofSM.sse - gofTCM.sse)/(gofTCM.sse/df2);
probF = fcdf(F, df1, df2, 'upper');
iFTest = probF<0.05;
dualfprintf(fid, ...
    ['Probability of observing by random chance the difference in sum-of-squares\n', ...
    'errors for the two-component model relative to that for the spherical model.\n\n']);
dualfprintf(fid,sprintf('P(F) = %.0f%%  ', probF*100));
if iFTest 
    dualfprintf(fid,'<<Two-component model has signficantly better fit>>\n\n')
else
    dualfprintf(fid,'<<Two-component model lacks signficantly better fit>>\n\n')
end
%... Assemble results for preferred fit
if iFTest
    %... Select results from two-component model
    coeff = coeffvalues(fitResultTCM);
    Beta = coeff(1);
    Rn = coeff(2);
    Rs = coeff(3);
    fn = coeff(4);
    residualZeta = zeta - fitResultTCM(r);
    rPred = linspace(0, max(r), 1000);
    zetaPred = fitResultTCM(rPred);
    rhoPred = tanh(zetaPred);
else
    %... Select results from spherical model
    coeff = coeffvalues(fitResultSM);
    DRho = coeff(1);
    Rs = coeff(2);
    Beta = 0;
    Rn = Rs;
    fn = 0;
    residualZeta = zeta - fitResultSM(r);
    rPred = linspace(0, max(r), 1000);
    zetaPred = fitResultSM(rPred);
    rhoPred = tanh(zetaPred);
end

dualfprintf(fid,'Best-Fit Grain Radius Results:\n');
dualfprintf(fid,sprintf('\t%3.2f (pixels)\n',Rs));
dualfprintf(fid,sprintf('\t%3.7f (mm)\n',Rs/Sample.pxPerMm));

%... Decompose into self and neighbor components
rPred = [0,linspace(min(r),max(r),1000)];
Sample.TCM.rPred = rPred;
Sample.TCM.rhoSC = (rPred<=2*Rs).*(1 - (3/4)*rPred/Rs + (1/16)*(rPred/Rs).^3);
Sample.TCM.rhoNC = fn*exp(-(rPred/Rn).^Beta).*sin(pi*rPred/Rn);


Sample.TCM.fit = fitResultTCM;
Sample.TCM.gof = gofTCM;
Sample.Rs = Rs;
Sample.Rs_micron = 1000*(Rs/Sample.pxPerMm);

end