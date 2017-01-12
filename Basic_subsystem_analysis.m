function [obj,y,c_ineq] = Basic_subsystem_analysis(subproblem_index,x_DV,PB)


% Constants
LAMBDA = PB.UserData.LAMBDA; % Default value is 0.

% default output values
obj = 0;
y = [];
c_ineq = [];

switch subproblem_index
    case 1
        % Subproblem 1
        [u,v,b] = get_variable(x_DV,PB,'u_1','v','b_1');  
        %fa = u.^1.5+1/(v+LAMBDA)+sqrt(b);
        a = log(u)+log(b)+log(v);
        fa = u+b+v+a;
        
        y = a;
        obj = fa;
    case 2
        % Subproblem 2
        [u,w,a] = get_variable(x_DV,PB,'u_2','w','a_2');            
        %b = log(1/LAMBDA+w^2-5*sqrt(w))+1/(u+LAMBDA)+1/(a+LAMBDA);
        %b = log(1/LAMBDA+w^2-sqrt(w))+1/(u+LAMBDA)+1/(a+LAMBDA);
        
        b = 1/(u+LAMBDA)+1/(w+LAMBDA)+1/(a+LAMBDA);
        y = b;
        c_ineq = w+b-10;
        obj = 0;

    otherwise
        error('unrecognized subproblem index')
end
