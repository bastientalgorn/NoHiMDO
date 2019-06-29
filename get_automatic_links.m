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
function L = get_automatic_links(i1,PB)
% Here, we automatically match 2 variables names if, when removing the trailing digits, the 
% 2 strings are identical and not empty.
% Example: "a_1" and "a_2" both give "a_" when removing trailing digits. That's a match.

% Links:
L = cell(0,0);

% Name of variable i1, without trailing spaces
s1 = remove_trailing_digits( PB.var_names{i1} );

% If s1 is empty, 
if isempty(s1)
    error(['Variable ' PB.var_names{i1} ' is empty after removing trailing digits']);
end

% Loop over all variables except i1
disp(['Variable ' PB.var_names{i1} ' is automatically linked with:']);
k = 1;
for i=1:PB.NV
    if (i==i1)
        continue;
    end
    si = remove_trailing_digits( PB.var_names{i} );
    if strcmp(s1,si)
        % Store link
        L{k} = i;
        k = k+1;
        % Display name
        disp(['    ' PB.var_names{i}]); 
    end 
end

if isempty(L)
    disp(['    Nothing...']);
end

% Display auto links:
function s = remove_trailing_digits(s)
while length(s) && isstrprop(s(end),'digit')
    s(end) = [];
end
