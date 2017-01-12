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
function v = get_value(pb_name,NO,NI,gamma,beta,optim_method,field)

dir = './matfiles/';

% Name of the file where to fetch the matrix M
mat_file = [dir pb_name '_' num2str(NO) '_' num2str(NI) '_' num2str(gamma) '_' num2str(beta) '_' optim_method '.mat'];

% FYI:
%     mat = [mean(tab_inc') ,  
%             std(tab_inc') ,
%             min(tab_inc') , 
%             max(tab_inc') , 
%        quantile(tab_inc',0.10) ,              
%        quantile(tab_inc',0.25) ,
%        quantile(tab_inc',0.50) ,         
%        quantile(tab_inc',0.75) ,
%        quantile(tab_inc',0.90) ,         
%            mean(tab_obj') ,  
%             std(tab_obj') , 
%             min(tab_obj') , 
%             max(tab_obj') , 
%        quantile(tab_obj',0.10) ,                
%        quantile(tab_obj',0.25) ,
%        quantile(tab_obj',0.50) ,           
%        quantile(tab_obj',0.75) ,
%        quantile(tab_obj',0.90) ];
%    mat = mat';

if findstr('inc',field)
    j = 0;
    div = 1;
elseif findstr('obj',field)
    j = 9;
    div = get_f_star(pb_name);
    
end

if findstr('mean',field)
    j = j+1;
elseif findstr('std',field)
    j = j+2;
elseif findstr('min',field)
    j = j+3;
elseif findstr('max',field)
    j = j+4;
elseif findstr('q0.10',field)
    j = j+5;
elseif findstr('q0.25',field)
    j = j+6;
elseif findstr('q0.50',field)
    j = j+7;
elseif findstr('median',field)
    j = j+7;
elseif findstr('q0.75',field)
    j = j+8;
elseif findstr('q0.90',field)
    j = j+9;
end

if findstr('curve',field)
    i = (1:NO);
else
    i = NO;
end



if exist(mat_file,'file')
    % Read the mat file
    load(mat_file)
    % Return the value.
    v = M(i,j)/div;
else
    % If no data, just produce a NaN matrix
    disp(mat_file);
    disp('File not found');
    v = nan(length(i),length(j));
end


