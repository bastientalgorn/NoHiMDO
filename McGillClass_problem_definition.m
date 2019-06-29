function PB = McGillClass_problem_definition(VersionNumber)



% There are 3 versions of the McGillClass problem.
% To allow for only one problem definition file and one
% analysis file, the version number will be store in the
% structure PB, in "UserData".
if nargin == 0
    VersionNumber = 3;
end
PB.UserData.VersionNumber = VersionNumber;

% The constant LAMBDA is also stored in PB
PB.UserData.LAMBDA = 0.1;


lb = 0;
ub = 10;

%            Name     SP      CV      links    dim    lb     ub
PB.var{1} = {'fu_1'   1    false      'fu_3'    1     lb     ub};
PB.var{2} = {'fa_1'   1    false      'auto'    1     lb     ub};
%----------------------------------------------------------
PB.var{3}  = {'a'     2    false         []     1     lb     ub};
PB.var{4}  = {'b'     2    false         []     1     lb     ub};
PB.var{5}  = {'c_2'   2    false         11     1     lb     ub};
PB.var{6}  = {'w_2'   2     true      'w_3'     1     lb     ub};
PB.var{7}  = {'fa_2'  2     true      'fa_1'    1     lb     ub};
%----------------------------------------------------------
PB.var{8}   = {'u'    3    false         []     1     lb     ub};
PB.var{9}   = {'v'    3    false         []     1     lb     ub};
PB.var{10}  = {'w_3'  3    false      'w_2'     1     lb     ub};
PB.var{11}  = {'c_3'  3     true      'c_2'     1     lb     ub};
PB.var{12}  = {'fu_3' 3     true     'fu_1'     1     lb     ub};
%----------------------------------------------------------
if VersionNumber==3
    % Add 2 variables to handle the inequality constraint in SP1
    PB.var{13}  = {'h'    1     true     'zero'     1    -100    +100};
    PB.var{14}  = {'zero' -1    true        'h'     1     0       0}; % Dummy
end


PB.index_main = 1;
PB.analysis_file = 'McGillClass_subsystem_analysis';

