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
function s = myformat(v,dummy)

if v == 0
    s = '0';
elseif v<0
    error('v < 0');
else
    L = floor(log10(v));
    if L>=0
        sign = '+';
    else
        if dummy
            sign = '--';
        else
            sign = '-';
        end
            
    end
    s = v/(10^L);
    L = num2str(abs(L));
    if length(L)==1
        L = ['0' L];
    end
    
    s = round(s*10);
    s = num2str(s);
    s = [s(1) '.' s(2)];
    s = [s 'e' sign L];
    
    e = (abs(v-str2num(s))/v);
    if e>0.1
        v
        s
        e
        error('error myformat');
    end

end

   
