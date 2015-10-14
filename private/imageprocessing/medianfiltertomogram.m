function [T] = medianfiltertomogram(Sample,T)
% this function will median filter a 2D image or 3D stack of images. The 3D
% stack of images is median filtered sequentially along each dimension

switch Sample.d
    case 3
        %... median filter each page
        h = waitbar(0,'Please Wait: Filtering Images'); 
        for i = 1:Sample.Nz;
            waitbar(i/Sample.Nz);    
            T(:,:,i) = medfilt2(T(:,:,i),Sample.medianFilter);
        end
        close(h);

        %... median filter along columns
        h = waitbar(0,'Please Wait: Filtering Images, orthogonal direction'); 
        for i = 1:Sample.Ny
            waitbar(i/Sample.Ny);    
            T(:,i,:) = medfilt2(squeeze(T(:,i,:)),Sample.medianFilter);
        end
        close(h);

        %... median filter along rows
        h = waitbar(0,'Please Wait: Filtering Images, final direction'); 
        for i = 1:Sample.Nx
            waitbar(i/Sample.Nx);    
            T(i,:,:) = medfilt2(squeeze(T(i,:,:)),Sample.medianFilter);
        end
        close(h);

    case 2
        %... median filter in 2D
        T = medfilt2(T,Sample.medianFilter);

    otherwise
        error('Specify 2 or 3 for Sample.d');
end
 

end