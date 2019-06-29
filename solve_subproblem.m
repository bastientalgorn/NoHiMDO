function subproblem_struct = solve_subproblem(PB,NoHi_options,nsp,x_target,x_master,delta_q,v,w,psize)

% Save rng state
if NoHi_options.save_subproblems
    subproblem_struct.rng_state = rng;
    subproblem_struct.nsp = nsp;
    subproblem_struct.x_target = x_target;
    subproblem_struct.x_master = x_master;
    subproblem_struct.delta_q = delta_q;
    subproblem_struct.v = v;
    subproblem_struct.w = w;
    subproblem_struct.psize = psize;
end

D = PB.D_indexes{nsp};
%C = PB.C_indexes{nsp};


% INIT
% Copy x_D from the x_master
x_D = x_master(D);
% The other components come from x_target
x = x_target;
x(D) = x_D;

%==================================================================
% Bounds
%==================================================================
lb = PB.lb(D);
ub = PB.ub(D);

%==================================================================
% Starting point(s) of the subproblem optimization
%==================================================================
list_starting_points = x_D;


%==================================================================
% Analysis functions
%==================================================================

%     function [fphi_c,y] = subproblem_analysis_overlay(x_D)
%         % Returns a vector fphi_c such that 
%         %    - 1st component is the objective value + penalty value
%         %    - next components are constraints values
%         % and a vector y containing the coupling value.
% 
%         % Compute the objective:
%         [f,y,constraints] = feval(PB.analysis_file,nsp,x_D,PB);
% 
%         % Affect Design and Coupling variables in x
%         x(D) = x_D;
%         if ~isempty(y)
%             x(C) = y;
%         end
% 
%         %Penality if coupling variables are outside their bounds.
%         if NoHi_options.constraints_cv
%             constraints = [ constraints , x(C)-ub(C) , lb(C)-x(C) ];
%         end
%         % Realistic objective
%         if NoHi_options.realistic_obj && nsp==PB.index_main
%             constraints = [constraints , f-PB.frealistic];
%         end
%         % Check for nan, inf and imaginary values
%         if any(isnan(x)) || any(imag(x))
%             f = +inf;
%         end
%         % Compute inconsistency
%         q = get_q(PB,x)+delta_q;
%         % Compute the consistency penalty:
%         % Multiply by Lagrange multipliers
%         phi = v.*q + w.*w.*q.*q;
%         % Sum on the relevant components of q
%         % (Note: it is possible to sum on all the components of pen
%         % but this perturbate the convergence when some values
%         % of pen are very large.)
%         phi = sum(phi(PB.Q_indexes{nsp}));
%         fphi_c = [f+phi constraints];
%         % Concatenate the non penalized objective at the end of y
%         % This allows to access the value via the "collect_y" option
%         % of simple mads.
%         y(end+1) = f;
%     end



%==================================================================
% Blackbox overlay allowing to slip the outputs 
% of the BB (objective, constraints, coupling variables)
%==================================================================

% fphi_c: vector such that:
%    - 1st component is the objective value + penalty value
%    - next components are constraints values
% y: coupling value returned by subsystem
saved_value_fphi_c = [];
saved_value_y = [];
x_D_last = [];

    function conditional_caller(x_D)
        if ~isequal(x_D,x_D_last)
            [saved_value_fphi_c,saved_value_y] = subproblem_analysis_overlay(PB,NoHi_options,nsp,x_D,x,v,w,delta_q);
            x_D_last = x_D;
        end
    end

    function fphi_out = function_fphi(x_D)
        conditional_caller(x_D);
        fphi_out = saved_value_fphi_c(1);
    end

    function [c_out,ceq_out] = function_c(x_D)
        conditional_caller(x_D);
        c_out = saved_value_fphi_c(2:end);
        ceq_out = [];
    end

    function y_out = function_y(x_D)
        conditional_caller(x_D);
        y_out = saved_value_y;
    end

%==================================================================
% Solvers
%==================================================================
switch NoHi_options.solver
    

    case 'mads'
        %==================================================================
        % MADS inputs
        %==================================================================
        f_handle = @(x_D) subproblem_analysis_overlay(PB,NoHi_options,nsp,x_D,x,v,w,delta_q);
        mads_options.budget=2*NoHi_options.NI*length(D);
        psize = psize(nsp);
        psize = max(psize,NoHi_options.tol*1000);
        mads_options.psize_init = min(1,psize); % Use the psize provided
        mads_options.check_cache = true;
        mads_options.tol = max(mads_options.psize_init/1000,NoHi_options.tol);
        mads_options.display = NoHi_options.solver_display;
        mads_options.collect_y = true;
        mads_options.opportunistic = false;
        mads_options.rich_direction= false;
        
        % Run optimization
        [x_D,~,mads_output] = simple_mads(list_starting_points,f_handle,lb,ub,mads_options);
        
%         disp('=======================================================');
%         disp(['SP ' num2str(nsp) ', nb success: ' num2str(mads_output.nb_success)]);
%         disp(['psize ratio: ' num2str(mads_options.psize_init/mads_output.psize) '  i.e.   ' ...
%             num2str(mads_options.psize_init) ...
%             ' --> ' num2str(mads_output.psize)]);
        
        % collect y_j from MADS output
        x_C = mads_output.ymin(1:end-1);
        if numel(x_C)~=numel(PB.C_indexes{nsp})
            nsp
            size(x_C)
            size(PB.C_indexes{nsp})
            whos
            error('x_C is not the right size');
        end
        % collect obj_j (the non penalised objective) from MADS output
        obj = mads_output.ymin(end);
        
        % compute the psize for the next iteration
        switch NoHi_options.psize
            case {'default','1',1}
                psize = 1;
            case 'success'
                psize = mads_output.psize_success;
            case 'max'
                psize = max( psize/16 , 2*mads_output.psize_max );
            case 'last'
                psize = mads_output.psize;
        end
        

    case {'sqp','interior-point','active-set','trust-region-reflective'}
        %==================================================================
        % fmincon
        %==================================================================

        fmincon_options = optimoptions('fmincon','Display','iter','Algorithm',NoHi_options.solver);
        
        % Run
        x_D = fmincon(@(x) function_fphi(x),x_D,[],[],[],[],lb,ub,@(x) function_c(x),fmincon_options);
        
        % collect y_j and obj_j (the non penalised objective)
        y = function_y(x_D);
        x_C = y(1:end-1);
        obj = y(end);
        
        % Dummu psize value
        psize = 0;
        
    otherwise
        error('NoHi_options.solver not recognized');
        
end


clear saved_value_fphi_c saved_value_y

%==================================================================
% Create subproblem_struct
%==================================================================
subproblem_struct.x_D = x_D;
subproblem_struct.x_C = x_C;
subproblem_struct.psize_out = psize;
subproblem_struct.obj = obj;

end