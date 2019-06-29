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

function [fphi_c,y] = subproblem_analysis_overlay(PB,NoHi_options,j,x_D,x,v,w,delta_q)

% Returns a vector fphi_c such that 
%    - 1st component is the objective value + penalty value
%    - next components are constraints values
% and a vector y containing the coupling value.

% Sets of indexes:
D = PB.D_indexes{j};
C = PB.C_indexes{j};
lb = PB.lb;
ub = PB.ub;

% Compute the objective:
% (This calls the analysis function defined by the user in "PB.analysis_file")
[f,y,constraints] = PB.analysis_handle(j,x_D,PB);

% Affect Design and Coupling variables in x
x(D) = x_D;
if ~isempty(y)
    x(C) = y;
end

%Penality if coupling variables are outside their bounds.
if NoHi_options.constraints_cv
    constraints = [ constraints , x(C)-ub(C) , lb(C)-x(C) ];
end

% Realistic objective
if NoHi_options.realistic_obj && j==PB.index_main
    constraints = [constraints , f-PB.frealistic];
end

% Check for nan, inf and imaginary values
if any(isnan(x)) || any(imag(x))
    f = +inf;
end

% Compute inconsistency
q = get_q(PB,x)+delta_q;
% Compute the consistency penalty:
% Multiply by Lagrange multipliers
phi = v.*q + w.*w.*q.*q;
% Sum on the relevant components of q
% (Note: it is possible to sum on all the components of pen
% but this perturbate the convergence when some values
% of pen are very large.)
phi = sum(phi(PB.Q_indexes{j}));

if isnan(phi) || imag(phi)
    phi = +inf;
end
fphi_c = [f+phi constraints];

% Concatenate the non penalized objective at the end of y
% This allows to access the value via the "collect_y" option
% of simple mads.
y(end+1) = f;
