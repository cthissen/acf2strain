function [] = acf2strain(Sample)
% This program calculates the best-fit deformation parameters from a 2-D
% image or a set of 3-D images. For given parameter values set in the
% Sample structure (see runsamples.m), the program will load the set of
% images, calculate the autocorrelation function (ACF) for the sample, and
% run a least squares fit for the deformation parameters that return the
% ACF to isotropy. The program also calculates estimated error parameters
% using the method of Hext (1963) and outputs both the results and a number
% of useful figures. For details of the method see: Thissen, C. J., &
% Brandon, M. T. (2015). An Autocorrelation Method for Three-Dimensional
% Strain Analysis. Journal of Structural Geology.
% 
% For comments, questions, or suggestions, please email cthissen@gmail.com
% or leave a comment under the issues tab at github.com/cthissen/acf2strain
% 
% Christopher J. Thissen, Yale University 
% Mark T. Brandon, Yale University
% Time-stamp: <Oct 09 2015>
% 
% Version history: 
% This version of the code is included as supplemental information associated with
% the following publication: Thissen, C. J., & Brandon, M. T. (2015). An 
% Autocorrelation Method for Three-Dimensional Strain Analysis. Journal of 
% Structural Geology. For the most up to date version, 
% visit github.com/cthissen/acf2strain or search the mathworks website.
% 
% v. 1.0 10/09/15: First release, Thissen and Brandon (Jour. of Structural 
% 				      Geology, 2015)
% v. 0.4 08/07/15: Updated file structure to use external functions
% v. 0.3 07/08/15: Added 2-component model for grain size estimate
% v. 0.2 05/19/15: removed weighting by lag
% v. 0.1 05/12/15: Added self-correlation grain size estimates
% v. 0.0 12/20/14: alpha version
% 
% REQUIREMENTS: This program can be RAM intensive due to the need to store many
% arrays of equivalent size to the three-dimensional tomogram. I have found that
% 16GB is sufficient for the examples included here. However, the code can be
% modified to run on systems with less RAM. Please contact the author for
% details.
% 
% This function requires the following toolboxes and additional codes: 
% 1. Curve Fitting Toolbox. 
% 2. Image Processing Toolbox. 
% 3. splinefit.m is included with the code, and is also available online at:
%    http://www.mathworks.com/matlabcentral/fileexchange/13812-splinefit
% 4. sort_nat is included with the code, and is also available online at:
%    http://www.mathworks.com/matlabcentral/fileexchange/10959-sort-nat--natural-
%     order-sort
% 5. export_fig is a useful suite of programs for saving publication-quality 
%    figures, and is included with the code. It is also available online at: 
%    http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% 6. coolwarm is included with the code. 
% 7. cmapscale is included with the code.
% 
% INSTALLATION: 
% No installation is required.
% 
% RUNNING The CODE: 
% We have included three examples: 1) the undeformed ooids, 2) the 3D
% phantom, 3) tomogram of sample 92810-3 from the Olympic Mountains. Each
% sample has a script with the prefix "run_" that can be used to run the
% strain analysis. For example, to run the strain analysis on the 3D
% phantom, open the "run_3Dphantom.m" file and push the "run" button, which
% looks like a green play button.
% 
% The required computation time for these examples will vary from computer
% to computer, but should be on the order of about 10 minutes.
% 
% OUTPUT: 
% The program will make a Results folder to store the output from the code,
% including best-fit deviatoric strain results and goodness of fit
% parameters. A number of useful images are plotted and saved, such as
% those found in Thissen and Brandon (2015). These can be found in the
% Results folder. If Flag.saveMatFile is set to true, the program will also
% save a.mat file in the Results folder. If present, this file will be
% loaded to avoid running the non-linear fitting routine again.

% Structure variables:
% Sample:           contains all information related to sample, including name,filename, etc
% Sample.flags:     contains all relevant information regarding code options, such as saving figures, movies, etc
% folder:           contains all relevant directory infomation, such as output folders and data folders
% StrainACF:        contains all relevant information regarding the code (version number etc)


% This code largely follows MATLAB programming conventions as specified here:
% (http://www.ee.columbia.edu/~marios/matlab/MatlabStyle1p5.pdf)

Folder = loadDirectoryInfo;

%% check for saved .mat file
fMat = [Folder.results,Sample.name,'.mat'];
chkExist = exist(fMat,'file');
switch chkExist
    case 2
        % mat file exists
        fprintf(1,'Previously saved mat file found, loading results. ');
        fprintf(1,'See file \n%s\nfor more information about sample results.\n',[Folder.results,Sample.name,'.txt']);
        %... store flag
        Flags = Sample.flag;
        load(fMat)
        Sample.flag = Flags;
        
        % reprint results
        fNameNew = [Folder.results,Sample.name,'_Reloaded.txt'];
        fprintf(1,'New results are saved in file:\n\t %s',fNameNew);
        Sample.fid = fopen(fNameNew,'w');
        printdeformationresults(Sample);
        
        % replot figures
        makeFigures(Sample,Folder);
        
        % remake movie
        if strcmp(Sample.flag.makeMovie,'true')
            makemovie(Sample)        
        end

        keyboard
        
        % exit function
        return
        
        
    case 0 
        % no matfile found
        fprintf(1,'No previously saved mat file found. Running program:\n');
end



%% continue
% prepare output file
Sample.fid = fopen([Folder.results,Sample.name,'.txt'],'w');
fid = Sample.fid;

% print information about program version
StrainAcf.version = 1.05;
dualfprintf(fid,sprintf('Sample Name: %s\n',Sample.name));
dualfprintf(fid,sprintf('Date: %s\n',date));
dualfprintf(fid,sprintf('StrainACF Version: %3.2f\n',StrainAcf.version));


%% Load and process tomogram

% Load tomogram
Sample = loadtomogram(Sample);
Sample.TOrig = int8(Sample.TOrig); % to save RAM, set original tomograph to 8 bit

% Crop tomogram 
Sample = croptomogram(Sample);
dualfprintf(Sample.fid,sprintf('Tomogram Cropped using Sample.cropParameters: %d %d %d %d\n',Sample.cropParameters));

% Ensure tomogram has odd dimensions so center of ACF grid is r=0
[Sample] = maketomogramdimensionsodd(Sample);
Sample.TcropOrig = Sample.T;

%% Make movie
% make movie of tomogram along sample z direction
cmap = gray(256);
Sample.cmapTomogram = cmapscale((Sample.T),cmap,0);
switch Sample.flag.makeMovie
    case 'false'
        % do nothing
    case 'true'
        makemovie(Sample)        
    otherwise
        error('Specify true or false for Flag.makeMovie');

end
close all
%% Print info about tomogram 
a = Sample.T;
a=whos('a');
MB = a.bytes/1048576;
dualfprintf(Sample.fid, sprintf('Image directory: \n\t %s \n',Sample.imagePath));
dualfprintf(Sample.fid, sprintf('Image Type: RBG = %i\n',   Sample.flag.rgb));
dualfprintf(Sample.fid, sprintf('Memory required for image = %2.2f MB\n',MB));
dualfprintf(Sample.fid, sprintf('Image Size (pixels) %d x %d x %d\n', Sample.Nx,Sample.Ny,Sample.Nz));
dualfprintf(Sample.fid, sprintf('Image resolution (pixels per mm) %3.2f\n', Sample.pxPerMm))
dualfprintf(fid,sprintf('Max pixel value r = : %3.2f\n',Sample.maxRadialLag));
dualfprintf(fid,sprintf('Data reduced to 1/2 sample volume by symmetry of ACF\n'));
dualfprintf(fid,sprintf('rho = 1 value has been removed from serial list of equations\n'));
dualfprintf(Sample.fid, sprintf('\n*******************\n'));
clear a

%% Center and Standardize Image

muT = mean(Sample.T(:));
dualfprintf(fid,sprintf('Calculating mean-centered ACF\n'));

% % also set "zingers" to mean
% Sample.T(Sample.T > 240) = muT;

sigmaT = sqrt(mean((Sample.T(:)-muT).^2));
Sample.Ts = (Sample.T-muT)/sigmaT;

%% Compute Autocorrelation Function

% get ACF function (5-7)
fprintf(1,'Calculating Autocorrelation Function\n');
Sample = computeautocorrelationfunction(Sample);

% Fisher transform data
Sample.zeta     = atanh(Sample.rho);
% Sample.zetaOrig = Sample.zeta;
Sample.zeta     = Sample.zeta(:);


%% Define grid and calculate weights
[Sample] = definesamplegrid(Sample);
[Sample] = reducedata(Sample);

% Precompute sample weights and degrees of freedom
[Sample.w, Sample.n]   = computeweights(Sample);
Sample.nu              = computedegreesoffreedom(Sample,Sample.n);        

%% Run Nonlinear Inversion 

% print info about isotropic ACF
dualfprintf(fid,sprintf('Fitting Function: SplineFit.m\n'));
dualfprintf(fid,sprintf('Cubic Polynomial.m\n'));        
dualfprintf(fid,sprintf('Number of Pieces: %d\n',(Sample.knots)));
       
% parameter order is defined as in Hext 1963 equation 2.1 (and 8.4)
%... E(1,1), E(2,2), E(3,3), E(2,3), E(3,1), E(2,1)
P0 = [0,0,0,0,0,0]; % initial guess for Henky tensor

% define options for non-linear fitting routine
options = optimset('TolFun',1e-6,'TolX',1e-6,'TolCon',1e-6,...
                   'LargeScale','on','Display','off','FinDiffType','central',...
                   'MaxFunEvals',1000);
% set bounds on fitting parameters               
upperBoundE =  1.5*ones(size(P0));
lowerBoundE = -1.5*ones(size(P0));

% define isochoric constraint function
switch Sample.d
    case 3
        % 3D tomogram
        isochoricFnc = @isochoricConstraint_3D;
         
    case 2
        % 2D image
        isochoricFnc = @isochoricConstraint_2D;
 
    otherwise
        error('Sample.d must be set to 2 or 3');
end


% run non-linear inversion
%... non-linear inversion subfunctions are setup to not modify anything in Sample structure
tic
[Sample.Pbest,Sample.Stilde,exitflag,output,lambda,grad,calcHessian] = ...
    fmincon(@(strain)chiSquared(strain,Sample),P0,[],[],[],[],lowerBoundE,upperBoundE,isochoricFnc,options);
timeFit=toc;

% Print some information about the fitting routine
dualfprintf(fid,'\n ********Best Fit Found**********\n');
dualfprintf(fid,sprintf('Time for best fit: %3.2f (sec)\n',timeFit));
dualfprintf(fid,sprintf('Function Evaluations: %d\n',output.funcCount));
dualfprintf(fid,sprintf('%s\n',output.message));


%% Calculate Grainsize, Best-fit, Residuals, and Covariance matrix

% Create best-fit solution from best-fit parameters
[Sample.Vhat,Sample.Ehat]   = straintotensor(Sample.Pbest);
[Sample.r0]                 = undeformingfunction(Sample,Sample.Ehat);   % get istropic coordinates using bestfit results
[Sample.zeta0,Sample.chat]  = computezeta0(Sample.r0,Sample);            % get best-fit isotropic function
[Sample]                    = estimategrainsize(Sample);

% keyboard
%... use 2-component model for best-fit results
% Sample.zeta0 = Sample.TCM.fit(Sample.r0);

% recompute weights and degrees of freedom (because they change in the case of cutoffLag weightingScheme).
[Sample.w,Sample.n] = computeweights(Sample);
[Sample.nu]         = computedegreesoffreedom(Sample,Sample.n);

% compute estimate of standard deviation
Sample.sigma0 = sqrt(Sample.Stilde);

% compute covariance matrix
Sample.H = Sample.nu/Sample.Stilde * calcHessian;
alpha = (Sample.H);
Sample.covEhat = inv(alpha);

switch Sample.d
    case 2
        dualfprintf(fid,'adjusting covariance matrix for 2D\n');
        Sample.covEhat(3:5,:) = 0;
        Sample.covEhat(:,3:5) = 0;
        Sample.covEhat;
    otherwise
        % no adjustment necessary
end

% calculate iSigma, the weight for each value 
Sample.iSigma = Sample.sigma0 ./ (Sample.w);

% calculate residuals
Sample.epsilonHat.raw                    = Sample.zeta - Sample.zeta0;
Sample.epsilonHat.unweightedStandardized = (Sample.zeta-Sample.zeta0)./Sample.sigma0;
Sample.epsilonHat.standardized           = Sample.w(:) .* (Sample.zeta-Sample.zeta0)./Sample.sigma0;

%% Assessment of Fit

% calculate R^2
Sample.R2 = 1-var(Sample.zeta-Sample.zeta0,1)./var(Sample.zeta(:),1);
% Sample.R2 = Sample.TCM.gof.rsquare;

% calculate durbin-watson statistic
[~,idx] = sort(Sample.r0,'ascend');
eHatAscending = Sample.epsilonHat.unweightedStandardized(idx); % this is the unweighted version
Sample.dw = sum(diff(eHatAscending).^2)/sum(eHatAscending(2:end).^2);

%% Hext Analysis: 

% Covariance and Confidence Intervals for Estimated Strain Parameters
stat = hextstatistics(Sample.Ehat,Sample.covEhat,Sample.nu,Sample.pValue);
Sample.defParams = henkytostretchstatistics(stat);


%... correct to positive plunge (NED reference frame)
for i = 1:3
    if Sample.defParams.vector{i}.plungeDegrees < 0
       Sample.defParams.vector{i}.plungeDegrees = -Sample.defParams.vector{i}.plungeDegrees;
       Sample.defParams.vector{i}.trendDegrees = Sample.defParams.vector{i}.trendDegrees + 180;
    end
end

% Print output to screen and file
printdeformationresults(Sample);
        
%% Save Mat File
switch Sample.flag.saveMatFile  
    case 'true'
        save([Folder.results,Sample.name,'.mat'],'-v7.3');
    case 'false'
        
    otherwise
        warning('Specify true or false for Flag.saveMatFile');
end
% keyboard

%% Make Figures
makeFigures(Sample,Folder)
 
%% Pause execution
% keyboard

end



