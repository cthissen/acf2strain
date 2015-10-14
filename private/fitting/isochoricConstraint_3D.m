function [c,ceq] = isochoricConstraint_3D(E)
    c = [];
    ceq = E(1) + E(2) + E(3); % require trace of Henky tensor is equal to 0
end