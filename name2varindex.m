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

function varindex = name2varindex(varname,PB)


try
    varindex = PB.name_to_index_struct.(varname);
catch
    % Use Levenshtein distance to find the closest var name
    disp(varname)
    PB
    d = +inf*ones(PB.NV,1)
    for i=1:PB.NV
        d(i) = levenshtein_dist(varname,PB.var_names{i});
    end
    dmin = min(d);
    disp(['Error: Cannot find any variable named "' varname]);
    
    s = ['Possible candidates: '];
    for i=1:PB.NV
        if d(i)==dmin
            s = [s '"' PB.var_names{i} '", '];
        end
    end
    disp( s(1:end-2) );

    error(['Variable ' varname ' could not be found.']);
end

%Levenshtein distance
function d=levenshtein_dist(s1,s2)
L1=numel(s1);
L2=numel(s2);
d=zeros([L2+1,L1+1]);
d(1,:)=0:L1;
d(:,1)=0:L2;
%Distance
for i=2:L2+1
    s2i=s2(i-1);
    for j=2:L1+1
        k = ~strcmp(s1(j-1),s2i);
        d(i,j) = min( [d(i-1,j-1)+k,d(i-1,j)+1,d(i,j-1)+1] );
    end
end
d = d(L2+1,L1+1);
