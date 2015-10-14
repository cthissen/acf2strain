function [Sample] = definesamplegrid(Sample)
% This function defines the grid for the image or tomogram. The is defined
% such that the center value is equal to (0,0,0) and the coordinates are
% integers.

% for the Sample structure, this function creates new fields x1,x2,x3, 
% which are vectors. The function also creates the Sample.r field, 
% which is radial distance


% rx, ry, and rz are defined using ndgrid so their matrix dimensions
%match the original image. Note however, that in the output ry and rx are
%flipped. This gives the intuitive results that rx increases as the columns
%increase (image right) and ry increases as the rows increase (image
%down). rz increases as the pages increase
rx = linspace(-(Sample.Nx-1)/2,(Sample.Nx-1)/2,Sample.Nx);
ry = linspace(-(Sample.Ny-1)/2,(Sample.Ny-1)/2,Sample.Ny);
rz = linspace(-(Sample.Nz-1)/2,(Sample.Nz-1)/2,Sample.Nz);
[ry,rx,rz] = ndgrid(rx,ry,rz);

Sample.rx = rx;
Sample.ry = ry;
Sample.rz = rz;


% get radial distance for each point
Sample.r = sqrt(Sample.rx.^2 + Sample.ry.^2 + Sample.rz.^2);


end

