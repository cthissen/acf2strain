% Run 3D phantom example
close all; clear all; clc

Sample.name            = 'Sx2.0_Sy1.0_Sz0.5_alpha30.0_10.02pxPerR';
Sample.sampleTitle     = 'Phantom Aggregrate'; % name printed on figures
Sample.imageFileExt    = '*.tif'; % extension for microtomography image files
Sample.medianFilter    = [];   % size of median filter
Sample.dimensions      = 3;    % number of sample dimensions
Sample.cropParameters  = [];   % parameters that define how to crop each slice
Sample.slices          = [1,209]; % min and max vertical slices to include
Sample.maxRadialLag    = 40;   % maximum radail lag to include in non-linear fit for deformation parameters
Sample.knots           = 100;  % number of pieces in spline fit for undeforming function
Sample.pxPerMm         = 100;  % number of pixels per millimeter
Sample.pValue          = 0.95; % confidence interval (two tailed)
 
Sample.imagePath    = ['Data',filesep,Sample.name,filesep]; % location of image files to read in

%... set preferences
Flag.medianFilter       = 'false';  % median filter removes "shot" noise from tomography
Flag.crop               = 'fullImage'; % automatic, fullImage, or define
Flag.plotFigures        = 'true';  % see plotfigures.m
Flag.makeMovie          = 'true';
Flag.saveFigures        = 'true'; 
Flag.saveMatFile        = 'true';


Sample.flag = Flag; % set as field in Sample

% run acf2strain
acf2strain(Sample)
    