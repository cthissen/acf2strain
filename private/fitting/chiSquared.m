function [Stilde] = chiSquared(strain,Sample)
% Given a candidate solution for the deformation parameters (store in the 
% variable, strain) this function calculates the objective function 
    
% parse strain vector into stretch and Hencky Tensors
[~,Etilde] = straintotensor(strain); 

[r0] = undeformingfunction(Sample,Etilde);

% calculate best-fit isotropic ACF
[zeta0Tilde] = computezeta0(r0,Sample);

 
% use precomputed weights
w = Sample.w;

% use precomputed degrees of freedom
nu = Sample.nu;


% objective function
Stilde = (1/nu) * sum(w(:) .* (Sample.zeta(:) - zeta0Tilde(:)).^2);


% sanity check
if isinf(Stilde)
    warning('isinf(Stilde)');    
    keyboard
elseif isnan(Stilde)
    warning('isnan(X2)');
    keyboard
elseif Stilde<0
    warning('Stilde<0');
    keyboard
elseif sum(imag(Stilde)>0)
    warning('imag(Stilde)')
    keyboard
end


% Write intermediate output to screen
fprintf(1,'Candidate Solution S: %3.8f \t E: % 3.8f % 3.8f % 3.8f % 3.8f % 3.8f % 3.8f \n',Stilde,strain);

end