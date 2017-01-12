function PB = Golinsky_problem_definition

f1min = 722;
f1max = 5408;
f2min = 184;
f2max = 506;
f3min = 942;
f3max = 1369;

%              Name    j     CV       tr      lb     ub
PB.var{1}  = {'x1_1'   1   false      2       2.6    3.6};
PB.var{2}  = {'x1_2'   2   false      [1,3]   2.6    3.6};
PB.var{3}  = {'x1_3'   3   false      2       2.6    3.6};

PB.var{4}  = {'x2_1'   1   false      5       0.7    0.8};
PB.var{5}  = {'x2_2'   2   false      [4,6]   0.7    0.8};
PB.var{6}  = {'x2_3'   3   false      5       0.7    0.8};

PB.var{7}  = {'x3_1'   1   false      8        17     28};
PB.var{8}  = {'x3_2'   2   false      [7,9]    17     28};
PB.var{9}  = {'x3_3'   3   false      8        17     28};

% SP1:
PB.var{10} = {'f1_1'   1   true       17    f1min  f1max};

% SP2:
PB.var{11} = {'f2_1'   2   true       18    f2min  f2max};
PB.var{12} = {'x4'     2   false      []      7.3    8.3};
PB.var{13} = {'x6'     2   false      []      2.9    3.9};

% SP3
PB.var{14} = {'f3_3'   3   true       19    f3min  f3max};
PB.var{15} = {'x5'     3   false      []      7.3    8.3};
PB.var{16} = {'x7'     3   false      []      5.0    5.5};

% SP4
PB.var{17} = {'f1_4'   4   false      10    f1min  f1max};
PB.var{18} = {'f2_4'   4   false      11    f2min  f2max};
PB.var{19} = {'f3_4'   4   false      14    f3min  f3max};


PB.index_main = 4;
PB.analysis_file = 'Golinsky_subsystem_analysis';


PB.fstar = 2994.47;
PB.frealistic = 2900;



