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

close all
clear all
disp('======= Solving the 3 problems for McGill Class ===============');

% NiHiMDO Parameters
NoHi_options.display = true;
NoHi_options.cache = false;
NoHi_options.w_scheme = 'median';
NoHi_options.beta = 1.5;
NoHi_options.NI = 30;
NoHi_options.NO = 30;
NoHi_options.save_subproblems = true;

PB = McGillClass_problem_definition(1);
output = NoHiSolver(PB,NoHi_options);
output

PB = McGillClass_problem_definition(2);
output = NoHiSolver(PB,NoHi_options);
output

PB = McGillClass_problem_definition(3);
output = NoHiSolver(PB,NoHi_options);
output
