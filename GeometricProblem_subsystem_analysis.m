function [obj,y,c_ineq] = GeometricProblem_subsystem_analysis(subproblem_index,x_DV)

% default values
obj = 0;
y = [];
c_ineq = [];


switch subproblem_index
    case 1
        % ATC - Subproblem 1
        [z3,z4,z5,z6,z7] = split_vector(x_DV);
        z1 = sqrt(z3^2+z4^-2+z5^2);
        z2 = sqrt(z5^2+z6^2 +z7^2);
        obj = z1^2+z2^2;
        c_ineq(1) = z3^-2 + z4^+2 - z5^+2;
        c_ineq(2) = z5^+2 + z6^-2 - z7^+2;
    case 2
        % ATC - Subproblem 2
        [z8,z9,z10,z11] = split_vector(x_DV);
        z3 = sqrt(z8^2+z9^-2+z10^-2+z11^2);
        y = z3;
        obj = 0;
        c_ineq(1) = z8^+2 + z9^+2  - z11^+2;
        c_ineq(2) = z8^-2 + z10^+2 - z11^+2;
    case 3
        % ATC - Subproblem 3
        [z11,z12,z13,z14] = split_vector(x_DV);
        z6 = sqrt(z11^2+z12^2+z13^2+z14^2);
        y = z6;
        obj = 0;
        c_ineq(1) = z11^+2 + z12^-2 - z13^+2;
        c_ineq(2) = z11^+2 + z12^+2 - z14^+2;
    otherwise
        error('unrecognized subproblem index')
end

