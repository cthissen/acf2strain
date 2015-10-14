% run sample 92810-3

Sample.name           = '92810-3';
Sample.sampleTitle    = Sample.name; % name printed on figures
Sample.imageFileExt   = '*.tif'; % extension for microtomography image files
Sample.medianFilter   = [3,3];   % size of median filter
Sample.dimensions     = 3;       % number of sample dimensions
Sample.cropParameters = [314.51  390.51  352.98  210.98];  % parameters that define how to crop each slice
Sample.slices         = [1,519]; % min and max vertical slices to include
Sample.maxRadialLag   = 67;      % maximum radail lag to include in non-linear fit for deformation parameters
Sample.knots          = 100;     % number of pieces in spline fit for undeforming function
Sample.pxPerMm        = 1000/7.54; % pixels per mm (see setup file)
Sample.pValue         = 0.95; % confidence interval (two tailed)


Sample.imagePath    = ['Data',filesep,Sample.name,filesep]; % location of image files to read in

%... set preferences
Flag.medianFilter       = 'true';  % median filter removes "shot" noise from tomography
Flag.crop               = 'automatic'; % automatic, fullImage, or define
Flag.plotFigures        = 'true';  % see plotfigures.m
Flag.makeMovie          = 'true';
Flag.saveFigures        = 'true'; 
Flag.saveMatFile        = 'true';

Sample.flag = Flag; % set as field in Sample

% run acf2strain
acf2strain(Sample)
    
