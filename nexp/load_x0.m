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

function x0 = load_x0(PB)

pb_name = PB.name;
!mkdir -p ./x0files
data_file = ['./x0files/x0_' pb_name '.mat'];


if ~exist(data_file,'file')
    NVAR = PB.NVAR;
    switch pb_name
        
        case 'biquad'
            NEXP = 40;
            x0 = randn(NEXP,NVAR);

        case 'geo'
            NEXP = 40;
            x0 = abs(randn(NEXP,NVAR));
            
        otherwise
            NEXP = 40;
            lb = PB.lb(:)';
            ub = PB.ub(:)';
            x0 = ones(NEXP,1)*lb + rand(NEXP,NVAR).*(ones(NEXP,1)*(ub-lb));

    end
    
    
    
    for k=1:NEXP
        for i=1:NVAR
            xi = x0(k,i);
            % What are the TR for this variable?
            J_tr = PB.tr{i};
            for j=J_tr
                xj = x0(k,j);

                x0(k,i) = 0.5*(xi+xj);
                x0(k,j) = x0(k,i);
            end
        end
    end
    
    
    save(data_file,'x0');
end



load(data_file);
