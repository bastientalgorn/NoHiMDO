%-------------------------------------------------------------------------------------%
%  NoHiMDO                                                                            %
%                                                                                     %
%  A solver for Multi-Disciplinary Optimization, based on Non-Hierarchical Analytical %
%  Target Cascading                                                                   %
%  Version 3.0.0                                                                      %
%                                                                                     %
%  Copyright (C) 2012-2019  Bastien Talgorn - McGill University, Montreal             %
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

function display_variables(PB,x)

NV = PB.NV;
NX = PB.NX;
NQ = PB.NQ;
NSP = PB.NSP;

% Display variable names of coupled variables
disp('Components of X:');
for i1=1:PB.NX
    s = ['x(' num2str(i1) ') = ' PB.x_names{i1}];
    if any(PB.link_matrix(i1,:))
        s = [s ' <---> { '];
        for i2=1:NX
            if PB.link_matrix(i1,i2)
                s = [s PB.x_names{i2} ', '];
            end
        end
        s = s(1:end-2);
        s = [s ' }'];
    end
    s = ['  ' s];
    disp(s);
end

disp('Variables:');
for i1=1:NV

    ximin = min(PB.var_index_to_x_indexes{i1});
    ximax = max(PB.var_index_to_x_indexes{i1});
    s = ['x(' num2str(ximin)];

    if ximin~=ximax
        s = [s ':' num2str(ximax)];
    end
    s = ['  ' s ') = ' PB.var_names{i1} ' = [ ' num2str(x(ximin:ximax)) ' ]'];
    disp(s);
end
   
% Display variable names of coupled variables
q = get_q(PB,x);
disp('Links:');
for k=1:PB.NQ
    i1 = PB.L(1,k);
    i2 = PB.L(2,k);
    s = ['q(' num2str(k) ') = ' PB.x_names{i1} ' - ' PB.x_names{i2} ' = ' num2str(q(k)) ' (SubPb. '];
    for j=1:NSP
        if ismember(k,PB.Q_indexes{j})
            s = [s num2str(j) ', '];
        end
    end
    s = s(1:end-2);
    s = ['  ' s '}'];
    disp(s)
end


