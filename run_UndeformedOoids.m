% run undeformed ooids example
close all; clear all; clc

Sample.name           = '2D_Crespi'; % folder name
Sample.sampleTitle    = 'Undeformed Ooids'; % name printed on figures
Sample.imageFileExt   = '*.tif'; % extension for microtomography image files
Sample.medianFilter   = [3,3];% size of median filter
Sample.dimensions     = 2; % number of sample dimensions
Sample.cropParameters =    1.0e+03 * [0.6215    0.1185    2.3750    2.0750]; % parameters that define how to crop each slice
Sample.slices         = [1,1];  % min and max vertical slices to include
Sample.maxRadialLag   = 1.5*241; % maximum radial lag to include in non-linear fit for deformation parameters
Sample.knots          = 100; % number of pieces in spline fit for undeforming function
Sample.pxPerMm        = 482.2; % number of pixels per millimieter
Sample.pValue         = 0.95; % confidence interval (two tailed)
Sample.scaleBarLength = 1; % scale bar length (in millimeters);

Sample.imagePath    = ['Data',filesep,Sample.name,filesep]; % location of image files to read in

%... set preferences
Flag.medianFilter       = 'true';  % median filter removes "shot" noise from tomography
Flag.crop               = 'define'; % automatic, fullImage, or define
Flag.plotFigures        = 'true';  % see plotfigures.m
Flag.makeMovie          = 'false';
Flag.saveFigures        = 'true'; 
Flag.saveMatFile        = 'true';

Sample.flag = Flag; % set as field in Sample


% run acf2strain
acf2strain(Sample)
    