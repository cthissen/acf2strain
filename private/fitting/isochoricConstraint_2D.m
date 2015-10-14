
function [c,ceq] = isochoricConstraint_2D(E)
    %... E(1,1), E(2,2), E(3,3), E(2,3), E(3,1), E(2,1)
    c = [];
    ceq(1) = E(1) + E(2); % require trace of Henky tensor = 0
    
    % require Z components to equal zero
    ceq(2) = E(3);
    ceq(3) = E(4);
    ceq(4) = E(5);

end