function [obj,y,c_ineq] = Impossible_subsystem_analysis(subproblem_index,x_DV,PB)

K = PB.UserData.K;

switch subproblem_index
    case 1
        % Subproblem 1
        [x,y] = get_variable(x_DV,PB,'x_1','y_1');  
        obj = norm([x+1 y+1]);
        c_ineq = x-y+K;
    case 2
        % Subproblem 2
        [x,y] = get_variable(x_DV,PB,'x_2','y_2');           
        obj = 0;
        c_ineq = norm([x-1 y])-1;
    otherwise
        error('unrecognized subproblem index')
end
y = [];