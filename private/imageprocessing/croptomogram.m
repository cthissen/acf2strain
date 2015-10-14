function [Sample] = croptomogram(Sample)
% This function crops a 2D or 3D stack of images. 
% If automatic crop is specified, function will not ask for new crop parameters.

% Replaces field T, Nx, Ny, and Nz with cropped version,
% defines Sample.cropParameters if not set earlier.


% get crop parameters
switch Sample.flag.crop
    case 'define'
        uiwait(msgbox('Draw a box around the location to crop all slices','Crop Tomogram or Image','modal'));
        
        % need to get crop parameters from user
        cropAtSlice = Sample.slices(1);
%         cropAtSlice = floor(Sample.Nz/2)+1;
        [~, Sample.cropParameters] = imcrop(Sample.T(:,:,cropAtSlice),gray(256));            

    otherwise
        % do nothing
end
            
        
switch Sample.flag.crop
    case 'fullImage'
        % do nothing

    otherwise
        % crop image
        hf=figure('visible','off'); %turns visibility of figure off 
        h = waitbar(0,'Please Wait: Cropping Images');

        %... initialize Tcropped
        Tcropped(:,:,1) = imcrop(Sample.T(:,:,1),Sample.cropParameters);
        Tcropped(:,:,2:Sample.Nz) = zeros([size(Tcropped(:,:,1)),Sample.Nz-1]);

        % crop image
        for i =1:Sample.Nz
             waitbar(i/Sample.Nz);
             Tcropped(:,:,i) = imcrop(Sample.T(:,:,i),Sample.cropParameters);
        end
        close(h);
        close(hf);   

        Sample.T = Tcropped; 
        
        %... only get slices we want
        Sample.T = Sample.T(:,:,Sample.slices(1):1:Sample.slices(2));
        
        [Sample.Nx,Sample.Ny,Sample.Nz] = size(Sample.T);
end
end
