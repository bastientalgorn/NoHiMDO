function PB = SBJ_problem_definition

addpath('./SBJ');

PB.UserData.altitude = 55000;
PB.UserData.Mach     = 1.4;



%             Name         SP     CV      link   dim   lb     ub   
PB.var{1}  = {'SFC_1'      1    false    'SFC_2'   1    1      4  };
PB.var{2}  = {'We_1'       1    false     'We_2'   1  100  30000  };
PB.var{3}  = {'LD_1'       1    false     'LD_3'   1  0.1     10  };
PB.var{4}  = {'Ws_1'       1    false     'Ws_4'   1 5000 100000  };
PB.var{5}  = {'Wf_1'       1    false         35   1 5000 100000  };
PB.var{6}  = {'Wt_1'       1     true         13   1 5000 100000  };

PB.var{7}  = {'D_2'        2    false         24   1 1000  70000  };
PB.var{8}  = {'T'          2    false         []   1  0.1      1  };
PB.var{9}  = {'We_2'       2     true          2   1  100  30000  };
PB.var{10} = {'SFC_2'      2     true          1   1    1      4  };
PB.var{11} = {'ESF_2'      2     true         12   1  0.5    1.5  };

PB.var{12} = {'ESF_3'      3    false         11   1  0.5    1.5  };
PB.var{13} = {'Wt_3'       3    false          6   1 5000 100000  };
PB.var{14} = {'theta_3'    3    false         36   1  0.2     50  };
PB.var{15} = {'tc_3'       3    false         28   1 0.01    0.1  };
PB.var{16} = {'ARw_3'      3    false         29   1  2.5      8  };
PB.var{17} = {'LAMBDAw_3'  3    false         30   1   40     70  };
PB.var{18} = {'Sref_3'     3    false         31   1  200    800  };
PB.var{19} = {'Sht_3'      3    false         32   1   50  148.9  };
PB.var{20} = {'ARht_3'     3    false         33   1  2.5    8.5  };
PB.var{21} = {'LAMBDAht'   3    false         []   1   40     70  };
PB.var{22} = {'Lw'         3    false         []   1 0.01    0.2  };
PB.var{23} = {'Lht'        3    false         []   1    1    3.5  };
PB.var{24} = {'D_3'        3     true          7   1 1000  70000  };
PB.var{25} = {'LD_3'       3     true          3   1  0.1     10  };
PB.var{26} = {'L_3'        3     true         27   1 5000 100000  };

PB.var{27} = {'L_4'        4    false         26   1 5000 100000  };
PB.var{28} = {'tc_4'       4    false         15   1 0.01    0.1  };
PB.var{29} = {'ARw_4'      4    false         16   1  2.5      8  };
PB.var{30} = {'LAMBDAw_4'  4    false         17   1   40     70  };
PB.var{31} = {'Sref_4'     4    false         18   1  200    800  };
PB.var{32} = {'Sht_4'      4    false         19   1   50  148.9  };
PB.var{33} = {'ARht_4'     4    false         20   1  2.5    8.5  };
PB.var{34} = {'Ws_4'       4     true          4   1 5000 100000  };
PB.var{35} = {'Wf_4'       4     true          5   1 5000 100000  };
PB.var{36} = {'theta_4'    4     true         14   1  0.2     50  };
PB.var{37} = {'lambda_4'   4    false         []   1  0.1    0.4  };
% Vectors for structural design
PB.var{38} = {'t'          4    false         []   9  0.1    4.0 };
PB.var{39} = {'ts'         4    false         []   9  0.1    9.0 };

PB.var{40} = {'zeros'     -1     true  'ineqtol'   1  0.0    0.0 };
PB.var{41} = {'ineqtol'    2     false  'zeros'   1  0.0    10.0 };



PB.index_main = 1;
PB.analysis_file = 'SBJ_subsystem_analysis';

PB.fstar = 33600;
PB.frealistic = 30000;
