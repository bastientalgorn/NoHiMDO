function PB = sac_problem_definition

lb = 0;
ub = 10;
ubW =1e+5;

%             Name   Sp      CV  link to   dim   lb     ub
PB.var{1} = {'Ww_1'   1    false   'Ww_2'    1    lb    ubW};
PB.var{2} = {'Wf_1'   1    false   'Wf_3'    1    lb    ubW};
PB.var{3} = {'x0_2'   2    false       6     2    lb     ub};
PB.var{4} = {'x2'     2    false      []     2    lb     ub};
PB.var{5} = {'Ww_2'   2     true       1     1    lb    ubW};
PB.var{6} = {'x0_3'   3    false       3     2    lb     ub};
PB.var{7} = {'x3'     3    false      []     2    lb     ub};
PB.var{8} = {'Wf_3'   3     true       2     1    lb    ubW};


PB.index_main = 1;


PB.analysis_file = 'sac_subsystem_analysis';

PB.fstar = 101027;
PB.frealistic = 80000;

