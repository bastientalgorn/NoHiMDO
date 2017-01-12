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

function [tab_inc,tab_obj,x] = NoHiSolver(PB,NoHi_options)


% Compute MDO struture.
if ~isfield(PB,'is_built') || ~PB.is_built
    PB = build_problem(PB);
end
display_variables(PB);


%==========================================================================
% Default options
%==========================================================================
default_NoHi_options.display = true;
default_NoHi_options.mads_display = false;
% Build an oracle based on the min value of
default_NoHi_options.oracle = false;
% Filter, after optimization, to go back to the least inconsistent design
default_NoHi_options.filter = false;
% Reduce the bounds of the optimization
default_NoHi_options.reduced_bounds = false;
% Updating scheme
default_NoHi_options.w_scheme = 'median';
% Build a cache for better initialization
default_NoHi_options.cache = true;
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
% Mads-Delta : special initialization of the poll size
default_NoHi_options.psize = 'success';
% Number of Inner/Outer loop iterations
default_NoHi_options.NI = 100;
default_NoHi_options.NO = 100;
% Hyper-parameters of the penalty update scheme
default_NoHi_options.beta = 2;
default_NoHi_options.gamma = 0.5;
% Initial value for the w vector
default_NoHi_options.w0 = 1;
% Initial design
default_NoHi_options.x0 = (PB.ub+PB.lb)/2;



%==========================================================================
% Use default options to complete options provided by user
%==========================================================================
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
% The filter method requires to build the CACHE
if NoHi_options.filter
    NoHi_options.cache = true;
end

%==========================================================================
% NHATC informations
%==========================================================================

% Number of outer loop iterations
NO = NoHi_options.NO;
% Number of sub-problem
NSP = PB.NSP;

% Starting point
x = NoHi_options.x0;
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
q = zeros(1,PB.NQ);

% Cache
global CACHE
if NoHi_options.cache
    if NoHi_options.display
        disp('Initialization of caches');
    end
    CACHE = cell(NSP,1);
    for nsp = 1:NSP
        CACHE{nsp} = [];
    end
end

% Result tabs
tab_obj = zeros(NO,1);
tab_inc = zeros(NO,1);
inc = 1;

%==========================================================================
% Set mads options
%==========================================================================

% Vector used to store the initial poll size
psize = ones(1,NSP);

for no = 1:NO

    if NoHi_options.display
        disp(['------------ Outer loop iteration ' num2str(no) '-----------']);
    end

    x_old = x;
    q_old = get_q(PB,x_old);
    
    %======================================================================
    % RUN SUB-PROBLEMS
    %======================================================================
    for nsp = 1:NSP
        D = PB.D_indexes{nsp};
        C = PB.C_indexes{nsp};
        V = PB.V_indexes{nsp};

        % INIT
        x_D = x(D);


        %==================================================================
        % Bounds
        %==================================================================
        lb = PB.lb(D);
        ub = PB.ub(D);
        if NoHi_options.reduced_bounds
            db = ub-lb;
            lb = max(lb,x_D-inc*db);
            ub = min(ub,x_D+inc*db);
        end

        %==================================================================
        % Starting point(s) of the subproblem optimization
        %==================================================================
        x0 = x_D;
        %---------------------------------
        % Oracle
        %---------------------------------
        if no>1 && NoHi_options.oracle
            x_oracle = x;
            q_opt = -v./(2*w.*w);
            for k=1:size(PB.L,2)
                i1 = PB.L(1,k);
                i2 = PB.L(2,k);
                if ismember(i1,D)
                    x_oracle(i1) = x_oracle(i2)+q_opt(k)*(PB.ub(i1)-PB.lb(i1));
                elseif ismember(i2,D)
                    x_oracle(i2) = x_oracle(i1)-q_opt(k)*(PB.ub(i1)-PB.lb(i1));
                end
            end
            x_oracle = x_oracle(D);
            % Snap to bounds
            x_oracle = min(x_oracle,PB.ub(D));
            x_oracle = max(x_oracle,PB.lb(D));
            % Create a set of candidate points between x and the oracle
            x0 = [];
            for z=linspace(0,0.01,5)
                x_try = (1-z)*x_D+z*x_oracle;
                x0 = [x0 ; x_try];
            end
        end
        %---------------------------------
        % Use cache to propose best candidate
        %---------------------------------
        if no>1 && NoHi_options.cache && size(CACHE{nsp},1)>0
            obj_pen = zeros(size(CACHE{nsp},1),1);
            x_try = x;
            for i=1:size(CACHE{nsp},1)
                xv = CACHE{nsp}(i,2:end);
                x_try(V) = xv;
                q_try = get_q(PB,x_try);
                obj_pen(i) = sum(v.*q_try+w.*w.*q_try.*q_try);
            end
            obj_pen = obj_pen+CACHE{nsp}(:,1);
            [dummy,i] = min(obj_pen);
            xv = CACHE{nsp}(i,2:end);
            x_try(V) = xv;
            x0 = [x_try(D);x0];
        end


        %==================================================================
        % Run optimization
        %==================================================================
        % FCN HANDLES
        f_handle = @(x_D) subproblem(PB,NoHi_options,nsp,x_D,x,v,w);
        mads_options.budget=2*NoHi_options.NI*length(D);
        mads_options.check_cache = true;       
        psize(nsp) = max(psize(nsp),NoHi_options.tol*1000);
        mads_options.psize_init = min(1,psize(nsp)); % Use the psize provided
        mads_options.tol = max(mads_options.psize_init/1000,NoHi_options.tol);    
        mads_options.opportunistic = false;
        mads_options.display = NoHi_options.mads_display;
                    
        [x_D,obj_pen,output] = simple_mads(x0,f_handle,lb,ub,mads_options);

        %==================================================================
        % compute the psize for the next iteration
        %==================================================================
        %disp(['psize = ' num2str(psize(nsp)) ' --> ' num2str(output.psize) ' ('  num2str(output.psize_success/output.psize) ' '  num2str(output.psize_max/output.psize) ')']);
        switch NoHi_options.psize
            case {'default','1',1}
                psize(nsp) = 1;
            case 'success'
                psize(nsp) = output.psize_success;
            case 'max'
                psize(nsp) = max( psize(nsp)/16 , 2*output.psize_max );
            case 'last'
                psize(nsp) = output.psize;
        end
        %psize(nsp) = max(mads_options.tol*10,psize(nsp));

        % Select best inc
        if NoHi_options.filter && no>1
            if NoHi_options.display
                disp('Filter');
            end
            %disp(['inc m : ' num2str(max(abs(get_q(PB,x_try)));
            inc = zeros(size(CACHE{nsp},1),1);
            x_try = x;
            for i=1:size(CACHE{nsp},1)
                xv = CACHE{nsp}(i,2:end);
                x_try(V) = xv;
                inc(i) = max(abs(get_q(PB,x_try)));
            end
            [inc_min,i] = min(inc);
            x_try(V) = CACHE{nsp}(i,2:end);
            x_D = x_try(D);
            if NoHi_options.display
                disp(['inc m : ' num2str(inc_min)]);
            end
        end

        % Perfrom ANALYSIS again for optimal design
        [obj_j,y] = feval(PB.analysis_file,nsp,x_D,PB);
        x(D) = x_D;
        x(C) = y;
    

        % If this is the main problem, then memorise the objective
        if nsp==PB.index_main
            obj = obj_j;
        end


    end


    %================================%
    % Compute the Inconsistencies    %
    %================================%
    q = get_q(PB,x);
    [inc,kmax] = max(abs(q));
    tab_obj(no:end) = obj;
    tab_inc(no:end) = inc;
    if NoHi_options.display
        % Absolute indexes of the two variables with the highest
        % inconsistency
        iabs1 = PB.L(1,kmax);
        iabs2 = PB.L(2,kmax);
        i1 = PB.x2v(iabs1);
        i2 = PB.x2v(iabs2);
        disp(['Highest inconsistency : ' PB.varnames{i1} ' : ' PB.varnames{i2} ]);
    end

    %======================================================================
    % Update v and w
    %======================================================================
    v = v+2.*w.*w.*q;
    switch NoHi_options.w_scheme
        case 'median'
            increase_w = (abs(q)>NoHi_options.gamma*(abs(q_old))) & (abs(q)>=median(abs(q)));
        case 'max'
            increase_w = (abs(q)>NoHi_options.gamma*(abs(q_old))) & (abs(q)>=max(abs(q)));
        case 'normal'
            increase_w = (abs(q)>NoHi_options.gamma*(abs(q_old)));
        otherwise
            disp(['NoHi_options.w_scheme = ' NoHi_options.w_scheme]);
            error('w_scheme not recognized');
    end
    w(increase_w) = w(increase_w)*NoHi_options.beta;


    %======================================================================
    % Display convergence
    %======================================================================
    if no==1
        dx = nan;
    else       
        dx = norm(x-x_old);
    end
    if true% NoHi_options.display
        %display_inconsistency(PB,x,x_old);
        Ls = 10;
        s = num2str(no); s(4) = ' ';
        s = [s ' || Inc: '];
        si = num2str(inc,3);
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

    %======================================================================
    % Stopping criterias
    %======================================================================
    if no>1 && abs(tab_inc(no))<NoHi_options.inc_stop 
        disp(['Stop: inc = ' num2str(inc,3) ' < ' num2str(NoHi_options.inc_stop,3) ]);
        break;
    end
    i = NoHi_options.noprogress_stop;
    if no>i+2 && log10(min(w))>6 && (tab_inc(no-i) <= min(tab_inc(no-i+1:no)))
        disp(['Stop: no progress after ' num2str(i) ' iterations.']);
        break;
    end

end



% Store value in PB.
PB.xfinal = x;
display_variables(PB);




