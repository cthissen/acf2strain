function [V,E]              = straintotensor(strain)
% reconstruct strain vector into stretch tensor and Hencky tensor
% parameter order is defined as in Hext 1963 equation 2.1 (and 8.4)
%... strain =  E(1,1), E(2,2), E(3,3), E(2,3), E(3,1), E(2,1)

    E(1,1) = strain(1);
    E(2,2) = strain(2);
    E(3,3) = strain(3); 
    
    E(2,3) = strain(4);
    E(3,1) = strain(5);
    E(2,1) = strain(6);
    
    E(1,2) = E(2,1);
    E(1,3) = E(3,1);
    E(3,2) = E(2,3);
    V = expm(E);
    
    
end