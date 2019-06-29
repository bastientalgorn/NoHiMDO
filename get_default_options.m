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

function default_NoHi_options = get_default_options

%==========================================================================
% Default options
%==========================================================================
default_NoHi_options.display = true;
% Updating scheme
default_NoHi_options.w_scheme = 'median';
% Initial value for the w vector
default_NoHi_options.constraints_cv = true;
% Realistic objective
default_NoHi_options.realistic_obj = false;
% Algo stops if this inconsistency is reached
default_NoHi_options.inc_stop = 1e-12;
% Stop criteria on psize
default_NoHi_options.tol = 1e-12;
% Algo stops if the inconsistency has not decreased for that many iter
default_NoHi_options.noprogress_stop = 100;
% Number of Inner/Outer loop iterations
default_NoHi_options.NI = 100;
default_NoHi_options.NO = 100;
% Hyper-parameters of the penalty update scheme
default_NoHi_options.beta = 2;
default_NoHi_options.gamma = 0.5;
% Initial value for the w vector
default_NoHi_options.w0 = 1;
% Initial design
default_NoHi_options.x0 = [];
% Nb of processors
default_NoHi_options.nb_proc = 1;
% Save detail of each subproblems
default_NoHi_options.save_subproblems = false;
% Solver
default_NoHi_options.solver = 'mads';
default_NoHi_options.solver_display = false;


% Mads options
%=========================================================

% Poll size update scheme 
% (How the initial poll size is computed from the final poll size of the previous subsystem optimization)
default_NoHi_options.psize = 'last';
