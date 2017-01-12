function [obj,y,c_ineq] = McGillClass_subsystem_analysis(subproblem_index,x_DV,PB)



% Constants
LAMBDA = PB.UserData.LAMBDA;
% ProblemNumber allow to know which of the 3 versions of this problem is
% used.
VersionNumber = PB.UserData.VersionNumber;
% Version 1 : Normal problem
% Version 2 : Same as 1 with one inequality constraint in SP2
% Version 3 : Same as 2 with one equality constriant in SP1


% default output values
obj = 0;
y = [];
c_ineq = [];

switch subproblem_index
    case 1
        % Subproblem 1
        [fu,fa] = get_variable(x_DV,PB,'fu_1','fa_1'); 
        obj = fu+fa;
        y = [];
        if VersionNumber==3
            h = fu-6*fa;
            y = h;
        end
    case 2
        % Subproblem 2
        [a,b,c] = get_variable(x_DV,PB,'a','b','c_2');  
        fa = a.^3+1/(b+LAMBDA)+LAMBDA*sqrt(c);
        w = log(1+fa+c);
        y = [w,fa];
        obj = 0;
    case 3
        % Subproblem 3
        [u,v,w] = get_variable(x_DV,PB,'u','v','w_3');        
        
        fu = sqrt(u+v)+1/(v+LAMBDA)+LAMBDA*w^2;
        c = log(11-w)+fu+1/(u+LAMBDA);
        y = [c,fu];
        if VersionNumber>=2
            c_ineq = v-3*u;
        end
        obj = 0;

    otherwise
        error('unrecognized subproblem index')
end
