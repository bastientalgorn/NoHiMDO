function PB = Impossible_problem_definition

PB.UserData.K = sqrt(2)-1+0.1;

lb = -3;
ub = +3;

%            Name     SP      CV      links    dim    lb     ub
%----------------------------------------------------------
PB.var{1}  = {'x_1'   1    false      'x_2'     1     lb     ub};
PB.var{2}  = {'y_1'   1    false      'y_2'     1     lb     ub};
%----------------------------------------------------------
PB.var{3}  = {'x_2'   2    false      'x_1'     1     lb     ub};
PB.var{4}  = {'y_2'   2    false      'y_1'     1     lb     ub};

% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 1;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'Impossible_subsystem_analysis';
PB.end_of_iter_file = 'Impossible_display';



