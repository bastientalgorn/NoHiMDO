function [obj,y,c_ineq] = Golinski_subsystem_analysis(subproblem_index,x_DV)

% default values
obj = 0;
y = [];
c_ineq = [];


switch subproblem_index
    case 1
        [x1,x2,x3] = split_vector(x_DV);
        % SP objective
        F1 = 0.7854*x1*x2^2*(3.3333*x3*x3 + 14.9335*x3 - 43.0934);
        y = F1;
        % Constraints
        g5 = 27/(x1*x2^2*x3) -1;
        g6 = 397.5/(x1*x2^2*x3^2) -1;
        g9 = x2*x3/40 -1;
        g10 = 5*x2/x1 -1;
        g11 = x1/(12*x2) -1;
        c_ineq = [g5,g6,g9,g10,g11];
    case 2
        [x1,x2,x3,x4,x6] = split_vector(x_DV);
        % SP objective
        F2 = -1.5079*x1*x6^2;
        F4 = 7.477*x6^3;
        F6 = 0.7854*x4*x6^2;
        y = F2+F4+F6;
        % Constraints
        g1 = sqrt( ((745*x4)/(x2*x3))^2 + 1.69e+7)/(110*x6^3) -1;
        g3 = (1.5*x6 + 1.9)/x4 -1;
        g7 = 1.93*x4^3/(x2*x3*x6^4) -1;
        c_ineq = [g1,g3,g7];
    case 3
        [x1,x2,x3,x5,x7] = split_vector(x_DV);
        % SP objective
        F3 = -1.5079*x1*x7^2;
        F5 = 7.477*x7^3;
        F7 = 0.7854*x5*x7^2;
        y = F3+F5+F7;
        % Constraints
        g2 = sqrt( ((745*x5)/(x2*x3))^2 + 1.575e+8)/(85*x7^3) -1;
        g4 = (1.1*x7 + 1.9)/x5 -1;
        g8 = 1.93*x5^3/(x2*x3*x7^4) -1;
        c_ineq = [g2,g4,g8];
    case 4
        [f1,f2,f3] = split_vector(x_DV);
        obj = f1+f2+f3;

    otherwise
        error('unrecognized subproblem index')
end
