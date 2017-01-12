function PB = GeometricProblem_problem_definition

lb = 1e-6;
ub = 1e+6;

%            Name     SP       CV       link     lb     ub
PB.var{1}  = {'z3_1'   1    false          2     lb     ub};
PB.var{2}  = {'z3_2'   2     true          1     lb     ub};
PB.var{3}  = {'z4'     1    false         []     lb     ub};
PB.var{4}  = {'z5'     1    false         []     lb     ub};
PB.var{5}  = {'z6_1'   1    false          6     lb     ub};
PB.var{6}  = {'z6_3'   3     true          5     lb     ub};
PB.var{7}  = {'z7'     1    false         []     lb     ub};
PB.var{8}  = {'z8'     2    false         []     lb     ub};
PB.var{9}  = {'z9'     2    false         []     lb     ub};
PB.var{10} = {'z10'    2    false         []     lb     ub};
PB.var{11} = {'z11_2'  2    false         12     lb     ub};
PB.var{12} = {'z11_3'  3    false         11     lb     ub};
PB.var{13} = {'z12'    3    false         []     lb     ub};
PB.var{14} = {'z13'    3    false         []     lb     ub};
PB.var{15} = {'z14'    3    false         []     lb     ub};

PB.index_main = 1;
PB.analysis_file = 'GeometricProblem_subsystem_analysis';

PB.fstar = 17.5887696038;
PB.frealistic = 15;

