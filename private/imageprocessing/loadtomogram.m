function [Sample] = loadtomogram(Sample)
% This function will load a 2D image or a series of 3D images in numerical
% order, as is the typical way of storing a tomogram.
% The function can currently handle 2D RGB or grayscale images, and a 3D
% stack of grayscale images. 

% adds fields Nx,Ny,Nz, T, and flag.rgb to Sample

% list all files with extension defined by imageFileExt
images = dir([Sample.imagePath,Sample.imageFileExt]);

%... number of images
N = numel(images);

%... parse file names
for i = 1:N
    imNames{i} = images(i).name;
end

%... sort file names to proper order
[~,index] = sort_nat(imNames);

% check for RGB image (NxNx3) or CMYK (NxNx4)
Sample.flag.rgb = 0;
[n,m,l] = size(imread([Sample.imagePath,images(1).name]));
if l == 3 
    Sample.flag.rgb = 1;
elseif l == 4
    error('no functionality for CMYK images at this time');
end



% load image according to image type:
h = waitbar(0,'Please Wait: Loading Images');
if N==1 && ~Sample.flag.rgb
    
    % 2D grayscale image
    imgFile = [Sample.imagePath,images.name];
    T = imread(imgFile);        


elseif (N~=1 && ~Sample.flag.rgb)
    
    % 3D stack of grayscale images
    %... initialize T
    imgFile = [Sample.imagePath,images(1).name];   
    tmp = imread(imgFile); 
    T = zeros(n,m,N);
    T(:,:,index(1)) = tmp; clear tmp;

    %... loop through each image file and load into T
    for j = 2:N
         i = index(j);
         waitbar(j/N);
         imgFile = [Sample.imagePath,images(i).name];
         T(:,:,j)=imread(imgFile);
    end

elseif N==1 && Sample.flag.rgb
    
    % 2D RGB image, convert to grayscale
    imgFile = [Sample.imagePath,images.name];
    tmp = imread(imgFile);
    T = rgb2gray(tmp);

elseif N~=1 && Sample.flag.rgb
    % stack of RGB images 
    error('no functionality for 3D RBG image stack at this time');        

else
    error('Error: imageLoad. Unspecified image file type');
end
close(h);


% get image size
[Sample.Nx,Sample.Ny,Sample.Nz] = size(T);


% define sample dimensionality
if Sample.Nz == 1
    Sample.d = 2;
else
    Sample.d = 3;
end


% median filter sample, if flag is set to true
switch Sample.flag.medianFilter
    case 'true'
        T= medianfiltertomogram(Sample,T);
        
    case 'false'
        % do not filter
        
    otherwise 
        error('Specify true or false for Flag.medianFilter');
end

Sample.T = double(T);
Sample.TOrig = Sample.T;

end
