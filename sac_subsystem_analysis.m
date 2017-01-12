function [obj,y,c_ineq] = sac_subsystem_analysis(subproblem_index,x_DV,PB)

% default values
obj = 0;
y = [];
c_ineq = [];


switch subproblem_index
    case 1
        % ATC - Subproblem 1
        [Ww,Wf] = get_variable(x_DV,PB,'Ww_1','Wf_1');   
        obj = 60000 + Ww + Wf;
    case 2
        % ATC - Subproblem 2
        [x0,x2] = get_variable(x_DV,PB,'x0_2','x2');
        Ww = 4000*(1+norm(x0-1)^2)*(1+norm(x2-1)^2);
        y = Ww;
        %obj = 0; % Already defined (beginning of this file)
    case 3
        % ATC - Subproblem 3
        [x0,x3] = get_variable(x_DV,PB,'x0_3','x3');
        xs = 10*(x0+x3);
        eggholder_handle = @(x) -(x(2)+47) * sin(sqrt(abs(x(2)+x(1)/2+47))) - x(1) * sin(sqrt(abs(x(1)-(x(2)+47))));
        EH = eggholder_handle(xs);
        omega = (1+norm(x0-2)^2)*(1+0.001*norm(x3-2)^2)*(1+1000*abs(EH));
        Dr = 0.025+0.004*log10(omega);
        Wf = 20000 + 380952*Dr + 9523809*Dr*Dr;
        y = Wf;
        %obj = 0; % Already defined (beginning of this file)
    otherwise
        error('unrecognized subproblem index')
end
