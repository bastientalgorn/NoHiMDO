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

function obj_pen = subproblem(PB,run_options,j,x_D,x,v,w)

% Sets of indexes:
D = PB.D_indexes{j};
C = PB.C_indexes{j};
V = PB.V_indexes{j};
lb = PB.lb;
ub = PB.ub;

% Compute the objective:
[obj,y,constraints] = feval(PB.analysis_file,j,x_D,PB);

% Affect Design and Coupling variables in x
x(D) = x_D;
if ~isempty(y)
    x(C) = y;
end

% Remember the evaluation
global CACHE
if run_options.cache
    if all(constraints<=0)
        CACHE{j} = [CACHE{j};[obj x(V)]];
    end
end

%Penality if coupling variables are outside their bounds.
if run_options.constraints_cv
    constraints = [ constraints , x(C)-ub(C) , lb(C)-x(C) ];
end

% Realistic objective
if run_options.realistic_obj && j==PB.index_main
    constraints = [constraints , obj-PB.frealistic];
end

% Check for nan, inf and imaginary values
if any(isnan(x)) || any(imag(x))
    obj = +inf;
end

% Compute the consistency penalty:
q = get_q(PB,x);
% Multiply by Lagrange multipliers
pen = v.*q + w.*w.*q.*q;
% Sum on the relevant components of q
% (Note: it is possible to sum on all the components of pen
% but this perturbate the convergence when some values
% of pen are very large.)
pen = sum(pen(PB.Q_indexes{j}));

% Returned value:
obj_pen = [obj+pen constraints];

