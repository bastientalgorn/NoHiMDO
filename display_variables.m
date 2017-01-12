%-------------------------------------------------------------------------------------%
%  NoHiMDO                                                                            %
%                                                                                     %
%  A solver for Multi-Disciplinary Optimization, based on Non-Hierarchical Analytical %
%  Target Cascading                                                                   %
%  Version 2.0.1                                                                      %
%                                                                                     %
%  Copyright (C) 2012-2016  Bastien Talgorn - McGill University, Montreal             %
%                                                                                     %
%  Author: Bastien Talgorn                                                            %
%  email: bastientalgorn@fastmail.com                                                 %
%                                                                                     %
%  This program is free software: you can redistribute it and/or modify it under the  %
%  terms of the GNU Lesser General Public License as published by the Free Software   %
%  Foundation, either version 3 of the License, or (at your option) any later         %
%  version.                                                                           %
%                                                                                     %
%  This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   %
%                                                                                     %
%  You should have received a copy of the GNU Lesser General Public License along     %
%  with this program. If not, see <http://www.gnu.org/licenses/>.                     %
%                                                                                     %
%  You can find information on NoHiMDO at https://github.com/bastientalgorn/NoHiMDO   %
%-------------------------------------------------------------------------------------%

function display_variables(PB)

disp('================ Summary of variables =========');
for i1 = 1:PB.NV
    
    % Variable index
    disp(['Variable #' num2str(i1)]);

    % Name
    disp(['  Name:        ' PB.varnames{i1}]);
    
    % Subproblem index
    dummy_str = '';
    if PB.dummy(i1)
        dummy_str = ' (dummy variable)';
    end
    disp(['  Subproblem:  ' num2str(PB.SP_index(i1)) dummy_str]);
    
    % Coupling variable?
    if PB.cv(i1)
        disp('  Coupl. var.: yes');
    else
        disp('  Coupl. var.: no');
    end
    
    % Links
    s = '  Linked to:  ';
    if isempty(PB.links{i1})
        s = [s ' nothing'];
    else
        for i2=PB.links{i1}
            s = [s ' ' PB.varnames{i2}];
        end
    end
    disp(s);
    
    % Dimension
    disp(['  Dimension:   ' num2str(PB.dim(i1))]);
    
    % Bounds
    disp(['  Lower bound: ' v2str(PB.lb(PB.X_indexes{i1}))]);
    disp(['  Upper bound: ' v2str(PB.ub(PB.X_indexes{i1}))]);
    
    % Final value (if available)
    if isfield(PB,'xfinal')
        disp(['  Solution:    ' v2str(PB.xfinal(PB.X_indexes{i1}))]);    
    end
end



function s = v2str(v)
if length(v)==1
    s = num2str(v,9);
elseif min(size(v))>1
    error('Not a vector');
else
    s = [ '[ ' num2str(v(:)',9) ' ]' ];
end
    
