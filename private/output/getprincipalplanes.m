function [Sample] = getprincipalplanes(Sample)
% this function creates planes normal to the best-fit principal strain axes

% get deformed coordinates
%... need to calculate these again because rx,ry,rz have been reduced
x1Vec = linspace(-(Sample.Nx-1)/2,(Sample.Nx-1)/2,Sample.Nx);
x2Vec = linspace(-(Sample.Ny-1)/2,(Sample.Ny-1)/2,Sample.Ny);
x3Vec = linspace(-(Sample.Nz-1)/2,(Sample.Nz-1)/2,Sample.Nz);
[x2,x1,x3] = ndgrid(x1Vec,x2Vec,x3Vec);   


Sample.T = Sample.T;


V = Sample.Vhat;
switch Sample.dimensions
    case 2
        % get undeformed coordinates    
        X1 = V(1,1)*x1 + V(1,2)*x2 + V(1,3)*x3;
        X2 = V(2,1)*x1 + V(2,2)*x2 + V(2,3)*x3;
%         X3 = V(3,1)*x1 + V(3,2)*x2 + V(3,3)*x3;
 
        % get undeformed acf and images
        SxNorm.plane.rho.undeformed = interp2(x1,x2,Sample.rho,X1,X2);
        SxNorm.plane.tomogram.undeformed = interp2(x1,x2,Sample.T,X1,X2);
 
    case 3

        % create multiples of principal directions that span Sample
        [vec,~] = sortEig(V,'descend');

        multiple = min([min(x1(:)),min(x2(:)),min(x3(:))]) : 1 : max([max(x1(:)),max(x2(:)),max(x3(:))]);
        vecParSx = [vec(1,1)*multiple; vec(2,1)*multiple; vec(3,1)*multiple];
        vecParSy = [vec(1,2)*multiple; vec(2,2)*multiple; vec(3,2)*multiple];
        vecParSz = [vec(1,3)*multiple; vec(2,3)*multiple; vec(3,3)*multiple];


        % make planes with axes parallel to principal directions
        %... planes are defined the principal strain direction normal to the plane
        h = waitbar(0,'Please wait: Making Principal Planes');
        for i = 1:numel(multiple)
            waitbar(i/numel(multiple));
            for j = 1:numel(multiple)
                SzNorm.x.deformed(i,j) = vecParSy(1,i) + vecParSx(1,j);
                SzNorm.y.deformed(i,j) = vecParSy(2,i) + vecParSx(2,j);
                SzNorm.z.deformed(i,j) = vecParSy(3,i) + vecParSx(3,j);

                SyNorm.x.deformed(i,j) = vecParSz(1,i) + vecParSx(1,j);
                SyNorm.y.deformed(i,j) = vecParSz(2,i) + vecParSx(2,j);
                SyNorm.z.deformed(i,j) = vecParSz(3,i) + vecParSx(3,j);

                SxNorm.x.deformed(i,j) = vecParSz(1,i) + vecParSy(1,j);
                SxNorm.y.deformed(i,j) = vecParSz(2,i) + vecParSy(2,j);
                SxNorm.z.deformed(i,j) = vecParSz(3,i) + vecParSy(3,j);

            end
        end
        close(h);

        % make undeformed planes
        SzNorm.x.undeformed = V(1,1)*SzNorm.x.deformed + V(1,2)*SzNorm.y.deformed + V(1,3)*SzNorm.z.deformed;
        SzNorm.y.undeformed = V(2,1)*SzNorm.x.deformed + V(2,2)*SzNorm.y.deformed + V(2,3)*SzNorm.z.deformed;
        SzNorm.z.undeformed = V(3,1)*SzNorm.x.deformed + V(3,2)*SzNorm.y.deformed + V(3,3)*SzNorm.z.deformed;

        SxNorm.x.undeformed = V(1,1)*SxNorm.x.deformed + V(1,2)*SxNorm.y.deformed + V(1,3)*SxNorm.z.deformed;
        SxNorm.y.undeformed = V(2,1)*SxNorm.x.deformed + V(2,2)*SxNorm.y.deformed + V(2,3)*SxNorm.z.deformed;
        SxNorm.z.undeformed = V(3,1)*SxNorm.x.deformed + V(3,2)*SxNorm.y.deformed + V(3,3)*SxNorm.z.deformed;

        SyNorm.x.undeformed = V(1,1)*SyNorm.x.deformed + V(1,2)*SyNorm.y.deformed + V(1,3)*SyNorm.z.deformed;
        SyNorm.y.undeformed = V(2,1)*SyNorm.x.deformed + V(2,2)*SyNorm.y.deformed + V(2,3)*SyNorm.z.deformed;
        SyNorm.z.undeformed = V(3,1)*SyNorm.x.deformed + V(3,2)*SyNorm.y.deformed + V(3,3)*SyNorm.z.deformed;

        % deformed Sample.Ts
        SxNorm.plane.tomogram.deformed = interp3(x1,x2,x3,Sample.T,SxNorm.x.deformed,SxNorm.y.deformed,SxNorm.z.deformed,'nearest');
        SyNorm.plane.tomogram.deformed = interp3(x1,x2,x3,Sample.T,SyNorm.x.deformed,SyNorm.y.deformed,SyNorm.z.deformed,'nearest');
        SzNorm.plane.tomogram.deformed = interp3(x1,x2,x3,Sample.T,SzNorm.x.deformed,SzNorm.y.deformed,SzNorm.z.deformed,'nearest');

        % undeformed Sample.Ts
        SxNorm.plane.tomogram.undeformed = interp3(x1,x2,x3,Sample.T,SxNorm.x.undeformed,SxNorm.y.undeformed,SxNorm.z.undeformed,'nearest');
        SyNorm.plane.tomogram.undeformed = interp3(x1,x2,x3,Sample.T,SyNorm.x.undeformed,SyNorm.y.undeformed,SyNorm.z.undeformed,'nearest');
        SzNorm.plane.tomogram.undeformed = interp3(x1,x2,x3,Sample.T,SzNorm.x.undeformed,SzNorm.y.undeformed,SzNorm.z.undeformed,'nearest');

        % store results for output
        Sample.SxNorm = SxNorm;
        Sample.SyNorm = SyNorm;
        Sample.SzNorm = SzNorm;
        
        Sample.vecParSx = vecParSx;
        Sample.vecParSy = vecParSy;
        Sample.vecParSz = vecParSz;

end

% store results for output

Sample.x1 = x1;
Sample.x2 = x2;
Sample.x3 = x3;



end