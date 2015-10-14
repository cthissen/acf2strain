function [misfit] = fminCRSforSelfCorrelation(params,Sample)
% this function returns the misfit for the self correlation for use in
% fmincrs to find the bounded best-fit parameters. 


Beta = params(1); 
   F = params(2);
Lambda=params(3);   
  RS = params(4); 
   Z = params(5);
r = Sample.r0(:);
 misfit = sum((((1-Z)*((r<2*RS).*(1 - (3/4)*r/RS + (1/16)*(r/RS).^3) ... 
     + (Lambda)*exp(-(r/(F*RS)).^Beta).*sin(2*pi*r/(F*RS)))+Z) ...
     - tanh(Sample.zeta(:))).^2);


fprintf('%05.3f \t%05.3f \t%05.3f \t%05.3f \t%05.3f \t%05.3f\n',misfit,params);

end