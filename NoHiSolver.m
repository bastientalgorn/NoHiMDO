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
%  Foundation, eitherBasi version 3 of the License, or (at your option) any later         %
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

function output = NoHiSolver(PB,NoHi_options)

%==========================================================================
% Build problem 
%==========================================================================
PB = build_problem(PB);


%==========================================================================
% Use default options to complete options provided by user
%==========================================================================
default_NoHi_options = get_default_options;

if ~exist('NoHi_options','var')
    NoHi_options = struct;
end
names = fieldnames(default_NoHi_options);
for i = 1:length(names)
    fieldstr = names{i};
    if ~isfield(NoHi_options,fieldstr)
        NoHi_options.(fieldstr) = default_NoHi_options.(fieldstr);
    end
end


%==========================================================================
% NHATC informations
%==========================================================================

% Number of outer loop iterations
NO = NoHi_options.NO;
% Number of sub-problem
NSP = PB.NSP;

% Starting point
if isempty(NoHi_options.x0)
    x = (PB.ub+PB.lb)/2;
else
    x = NoHi_options.x0;
end


% Check starting point
x = x(:)';
if length(x)<PB.NX
    warning('Correcting dimension of x');
    lx = length(x);
    x(lx+1:PB.NX) = default_NoHi_options.x0(lx+1:PB.NX);
end
if length(x)>PB.NX
    warning('Correcting dimension of x');
    x = x(1:PB.NX);
end

% Initialization of the penalty vectors
v = zeros(1,PB.NQ);
w = NoHi_options.w0*ones(1,PB.NQ);

% Result tabs
tab_obj = zeros(NO,1);
tab_inc = zeros(NO,1);

if NoHi_options.display
    display_variables(PB,x);
end

% Cell array to store the results of each subproblem
if NoHi_options.nb_proc ~= 1
    result = cell(1,NSP);
    result(1:NSP) = {struct};
else
    result = struct;
end

% Stopping criteria (Reason why the outer loop stopped).
stop = "";


% Vector used to store the initial poll size
psize = ones(1,NSP);

% Outer loop.
for no = 1:NO

    % Copy of values of x and q before iteration.
    x_old = x;

    %===========================================%
    %===========================================%
    %                                           %
    %          SOLVE SUB-PROBLEMS               %
    %                                           %
    %===========================================%
    %===========================================%
  
    if NoHi_options.nb_proc ~= 1
        %===================================
        % PARALLEL MODE
        %===================================
        % Compute target (averaged values of linked variables)
        x_target = x*(PB.xavg_matrix');
        
        % Compute correction delta_q for calculation of q
        delta_q = -v./(2*w.^2)-get_q(PB,x_target);
        
        % Solve subproblem in parallel
        parfor nsp=1:NSP
            result{nsp} = solve_subproblem(PB,NoHi_options,nsp,x_target,x,delta_q,v,w,psize);
        end
        % Update x
        for nsp = 1:NSP
            x(PB.D_indexes{nsp}) = result{nsp}.x_D;
            x(PB.C_indexes{nsp}) = result{nsp}.x_C;
            psize(nsp) = result{nsp}.psize_out;
        end
        % Store subproblem details in output.
        if NoHi_options.save_subproblems
            for nsp = 1:NSP
                output.subproblems(nsp,no) = result{nsp};
            end
        end
        obj = result{PB.index_main}.obj;
    else
        %====================
        % SERIAL MODE
        %=====================
        delta_q = 0;
        for nsp=1:NSP
            result = solve_subproblem(PB,NoHi_options,nsp,x,x,delta_q,v,w,psize);
            % Update x
            x(PB.D_indexes{nsp}) = result.x_D;
            x(PB.C_indexes{nsp}) = result.x_C;
            psize(nsp) = result.psize_out;
            % If this is the main problem, then memorise the objective
            if nsp==PB.index_main
                obj = result.obj;
            end
            if NoHi_options.save_subproblems
                output.subproblems(nsp,no) = result;
            end
        end
    end
    
    %===========================================%
    %===========================================%
    %                                           %
    %                DISPLAY                    %
    %                                           %
    %===========================================%
    %===========================================%

    %================================%
    % Compute the Inconsistencies    %
    %================================%
    q = get_q(PB,x);
    [qmax,kmax] = max(abs(q));
    tab_obj(no:end) = obj;
    tab_inc(no:end) = qmax;
    if NoHi_options.display
        % Absolute indexes of the two variables with the highest
        % inconsistency
        i1 = PB.L(1,kmax);
        i2 = PB.L(2,kmax);
        disp(['Highest inconsistency : ' PB.x_names{i1} ' vs ' PB.x_names{i2} ]);
    end

    %======================================================================
    % User's display
    %======================================================================
    if isfield(PB,'end_of_iter_file')
        eval([PB.end_of_iter_file '(PB,x,v,w);']);
    end


    %======================================================================
    % Display convergence
    %======================================================================
    dx = norm(x-x_old);
    if NoHi_options.display
        Ls = 10;
        s = num2str(no); s(4) = ' ';
        s = [s ' || qmax: '];
        si = num2str(qmax,3);
        si(Ls) = ' '; si = si(1:Ls); s = [s si];
        s = [s ' || Obj: '];
        si = num2str(obj,3);
        si(Ls) = ' '; si = si(1:Ls); s = [s si];
        s = [s ' || dx: '];
        si = num2str(dx,3);
        si(Ls) = ' '; si = si(1:Ls); s = [s si];
        s = [s ' || max(w): '];
        si = num2str(max(w),3);
        si(Ls) = ' '; si = si(1:Ls); s = [s si];
        disp(s);
    end

  
    %===========================================%
    %===========================================%
    %                                           %
    %                UPDATES                    %
    %                                           %
    %===========================================%
    %===========================================%


    %======================================================================
    % Update v and w
    %======================================================================
    v = v+2.*w.*w.*q;
    q_old = get_q(PB,x_old);
    q_stall = (abs(q)>NoHi_options.gamma*(abs(q_old)));
    switch NoHi_options.w_scheme
        case 'median'
            increase_w = 2 * (q_stall & (abs(q)>=median(abs(q))));
        case 'max'
            increase_w = length(q) * (q_stall & (abs(q)>=max(abs(q))));
        case 'normal'
            increase_w = q_stall;
        case 'rank'
            [~,rank] = sort(q);
            increase_w = 2 * q_stall .* rank/max(rank);
        otherwise
            disp(['NoHi_options.w_scheme = ' NoHi_options.w_scheme]);
            error('w_scheme not recognized');
    end
    w = w.*(NoHi_options.beta.^increase_w);

    %======================================================================
    % Stopping criterias
    %======================================================================
    if no>1 && abs(tab_inc(no))<NoHi_options.inc_stop 
        disp(['Stop: qmax = ' num2str(qmax,3) ' < ' num2str(NoHi_options.inc_stop,3) ]);
        stop = "Max inconsitency is below stopping threshold";
        break;
    end
    i = NoHi_options.noprogress_stop;
    if no>i+2 && log10(min(w))>6 && (tab_inc(no-i) <= min(tab_inc(no-i+1:no)))
        disp(['Stop: no progress after ' num2str(i) ' iterations.']);
        stop = "No progress after several iterations";
        break;
    end

end

if isempty(stop)
    stop = "Iteration budget exhausted";
end
    

% Store value in PB.
if NoHi_options.display
    %display_variables(PB,x);
end

% Store outputs
output.x = x;
output.obj = obj;
output.tab_inc = tab_inc;
output.tab_obj = tab_obj;
output.stop = stop;
output.q = q;
 

