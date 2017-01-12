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

function myboxplot(x,vmean,vmin,vmax,vq1,vq2,vq3,color,dir)

L = 0.4;

if nargin<=7
    color = 0.9*[1 1 1];
end
if nargin<=7
    dir = 1;
end

hold on;

vl = vq1;
vu = vq3;

if dir==1
    fill(x+L*[-1 +1 +1 -1 -1],[vl vl vu vu vl],color);
    plot(x*[1 1],[vmin vmax],'k');
    plot(x+[-1 1]*L/5,vmin*[1 1],'k');
    plot(x+[-1 1]*L/5,vmax*[1 1],'k');
    plot(x+[-L +L],vq2*[1 1],'k','linewidth',3);
    %plot(x,vmean,'ok');
else
    fill([vl vl vu vu vl],x+L*[-1 +1 +1 -1 -1],color);
    plot([vmin vmax],x*[1 1],'k');
    plot(vmin*[1 1],x+[-1 1]*L/5,'k');
    plot(vmax*[1 1],x+[-1 1]*L/5,'k');
    plot(vq2*[1 1],x+[-L +L],'k','linewidth',3);
    %plot(x,vmean,'ok');
end



