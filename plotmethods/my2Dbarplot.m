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
function my2Dbarplot(m)

[NI,NJ] = size(m);

text_size = 12;
mode = 'detached';
%  mode = 'grouped';
%  mode = 'stacked';

bar3(m,mode)
set(gca,'ZScale','log')

llim = get(gca,'zlim');
llim = llim(1);


h = get(gca,'Children');

for i = 1:length(h)
       ZData = get(h(i), 'ZData');
       ZData(ZData==0) = llim;
       set(h(i), 'ZData', ZData);
       set(h(i),'facecolor',(0.7+0.2*i/length(h))*[1 1 1]);
       set(h(i),'edgecolor',0.4*[1 1 1]);
end

%dm = max(max(m))-min(min(m));

zl = get(gca,'zlim');
rz = (zl(2)/zl(1))^0.03;

for i=1:NI
    for j=1:NJ
        text(j,i,max(m(i,j),llim)*rz,myformat(m(i,j),false),'horizontalalignment','center','fontsize',text_size)
    end
end
