function [r0] = undeformingfunction(Sample,E)
% given a candidate solution E, calculate r0, the magnitude of the initial
% lag vector

Binv = expm(-2*E);

r0 = sqrt(Binv(1,1)*Sample.rx.^2 ...
         +Binv(2,2)*Sample.ry.^2 ...
         +Binv(3,3)*Sample.rz.^2 ...
         +2*Binv(1,2)*Sample.rx.*Sample.ry ...
         +2*Binv(1,3)*Sample.rx.*Sample.rz ...
         +2*Binv(2,3)*Sample.ry.*Sample.rz);
end