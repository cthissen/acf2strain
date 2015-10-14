function [Sample] = maketomogramdimensionsodd(Sample)
% this function makes the tomogram dimensions odd (so the center value is r=0)

% replaces fields T, Nx,Ny,Nz in Sample

% find largest odd value for each dimension
oddDim = [Sample.Nx,Sample.Ny,Sample.Nz];
even = mod(oddDim,2)==0;
oddDim(even) = oddDim(even)-1;
        


% apply crop to make odd
switch Sample.d
    case 2
    Sample.T = Sample.T(1:oddDim(1),1:oddDim(2));
        
    case 3
    Sample.T = Sample.T(1:oddDim(1),1:oddDim(2),1:oddDim(3));
    
end

[Sample.Nx,Sample.Ny,Sample.Nz] = size(Sample.T);

end