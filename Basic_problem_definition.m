function PB = Basic_problem_definition


% The constant LAMBDA is also stored in PB
PB.UserData.LAMBDA = 0;


lb = 0;
ub = 10;

%            Name     SP      CV      links    dim    lb     ub
%----------------------------------------------------------
PB.var{1}  = {'u_1'   1    false      'u_2'     1     lb     ub};
PB.var{2}  = {'v'     1    false        []      1     lb     ub};
PB.var{3}  = {'a_1'   1     true      'a_2'     1     lb     ub};
PB.var{4}  = {'b_1'   1    false      'b_2'     1     lb     ub};
%----------------------------------------------------------
PB.var{5}   = {'u_2'  2    false      'u_1'     1     lb     ub};
PB.var{6}   = {'w'    2    false        []      1     lb     ub};
PB.var{7}   = {'a_2'  2    false      'a_1'     1     lb     ub};
PB.var{8}   = {'b_2'  2     true      'b_1'     1     lb     ub};


PB.index_main = 1;
PB.analysis_file = 'Basic_subsystem_analysis';

