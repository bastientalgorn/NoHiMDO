%-------------------------------------------------------------------------------------%
%  simple_mads                                                                        %
%                                                                                     %
%  A simple matlab version of the Mesh Adaptive Direct Search algorithm               %
%  for constrained derivative free optimization.                                      %
%  Version 2.0.3                                                                      %
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
%  You can find information on simple_mads at                                         %
%  https://github.com/bastientalgorn/simple_mads                                      %
%-------------------------------------------------------------------------------------%

function [xmin,fmin,output] = simple_mads(X0,bb_handle,lb,ub,options)

%-------------------------------
% Options
%-------------------------------
% Default options
default_options.budget        = +inf;
default_options.tol           = 1e-9;
default_options.psize_init    = 1.0;
default_options.display       = false;
default_options.opportunistic = true;
default_options.check_cache   = true; % Do not evaluate the same point twice
default_options.store_cache   = false; % Store and return the list of candidates evaluated.
default_options.collect_y     = false;
default_options.rich_direction= true;

% Construct the option structure from user options and default options
if ~exist('options','var')
    options = struct;
end
names = fieldnames(default_options);
for i = 1:length(names)
    fieldstr = names{i};
    if ~isfield(options,fieldstr)
        options.(fieldstr) = default_options.(fieldstr);
    end
end

%-------------------------------
% Dimension and verifications
%-------------------------------
DIM = size(X0,2);
if length(lb)~=DIM
    error('Wrong lb dimension');
end
if length(ub)~=DIM
    error('Wrong ub dimension');
end

%-------------------------------
% Initialization
%-------------------------------
iteration = 0;
bb_eval = 0;
nb_success = 0;
fmin = +inf;
hmin = +inf;

% If the blackbox has 2 return arguments and if collect_y is activated, then the 2nd argument is store in ytry.
% The value of ytry associated with the best candidate xmin is store in ymin and returned.
ymin = [];
ytry = [];

% Variable scaling
scaling = (ub-lb)/10.0;
scaling(isinf(scaling)) = 1.0;
scaling = diag(scaling);

% Poll size initialization
psize = options.psize_init;
psize_success = 0;
psize_max = 0;

% Hash function 
hv = [1+rand 1+rand 1+rand 1+rand (1+rand)*1e-6 1+rand];
hashfcn = @(x) prod( (mod(x,hv(1))+hv(2)).*(mod(x,hv(3))+hv(4)).*(1e+6*mod(x,hv(5))+hv(6)) );
hashtable = [];
% Cache table
cache = [];


%-------------------------------
% Start the optimization
%-------------------------------
while true

    if iteration==0
        % Evaluate starting points
        if (options.display)
            disp('Evaluation of the starting points');
        end
        POLL = X0;
    else
        %-------------------------------
        % Build polling directions
        %-------------------------------
        msize = min(psize^2,psize);
        rho = psize/msize;
        if options.rich_direction
            % Generate direction
            v = randn(DIM,1);
            % Normalize
            v = v/norm(v);
            % Build Householder matrix
            H = eye(DIM)-2*v*v';
        else
            H = eye(DIM);
        end
        % Normalization of each column
        H = H*diag(max(abs(H)).^-1);
        % Rounding (and transpose)
        H = msize*ceil(rho*H)';
        % Add the opposite directions and scale
        H = [H;-H]*scaling;
        % Build POLL / central point
        POLL = bsxfun(@plus,xmin,H);
        % Shuffle poll
        POLL = POLL(randperm(2*DIM),:);
        % Snap to bounds
        POLL = bsxfun(@min,POLL,ub);
        POLL = bsxfun(@max,POLL,lb);
    end

    %-------------------------------
    % Evaluate points of the poll
    %-------------------------------
    success = false;
    for i=1:size(POLL,1)

        % Get point
        xtry = POLL(i,:);

        % Check if xtry has already been evaluated
        if options.check_cache && length(hashtable)
            xtry_hash = hashfcn(xtry);
            if ismember(xtry_hash,hashtable)
                if options.display
                    disp('Cache hit');
                end
                continue;
            else
                hashtable(end+1) = xtry_hash;
            end
        end

        % Evaluation of the blackbox
        if options.collect_y   
            [bb_output,ytry] = bb_handle(xtry);
        else
            bb_output = bb_handle(xtry);
        end
        bb_eval = bb_eval+1;
        
        % Objective function
        ftry = bb_output(1);
        % Constraints (can be an empty vector)
        ctry = bb_output(2:end);
        % Aggregate constraint
        htry = sum(max(ctry,0).^2);
        if isnan(htry) || any(isnan(ctry))
            htry = +inf;
        end
        % Penalize the objective
        if isnan(ftry) || (htry>0)
            ftry = +inf;
        end

        % Add to the cache
        if options.store_cache
            cache(end+1,:) = xtry;
        end

        % Test for success
        if ( (hmin>0) && (htry<hmin) ) || ( (htry==0) && (ftry<fmin) )
            success = true;
            nb_success = nb_success+1;
            xmin = xtry;
            fmin = ftry;
            hmin = htry;
            ymin = ytry;
            if options.display
                disp(['Succes: ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
            end
            psize_success = psize;
            psize_max = max(psize,psize_max);
        end
        % Test for stopping criteria
        if bb_eval>=options.budget
            break; % Reached the total budget
        end
        if success && options.opportunistic && (iteration>1)
            break; % Quit the evaluation of the POLL
        end

    end

    %-------------------------------
    % Updates
    %-------------------------------
    if iteration>0
        if success
            psize = psize*2;
        else
            psize = psize/2;
        end
    end

    if options.display
        if iteration==0
            disp('End of the evaluation of the starting points');
        end
        disp(['iteration=' num2str(iteration) ...
              ' bb_eval=' num2str(bb_eval) ...
              ' psize=' num2str(psize,3) ...
              ' hmin=' num2str(hmin,3) ...
              ' fmin=' num2str(fmin,3)]);
    end

    if (abs(psize)<options.tol) || (bb_eval>=options.budget)
        break;
    end

    iteration = iteration+1;
end % end of the mads iterations

if options.display
    disp('end of simple_mads');
    disp(['Final objective value: ' num2str(fmin) ' (hmin = ' num2str(hmin) ')']);
end

% Build the output
output.xmin = xmin;
output.fmin = fmin;
output.hmin = hmin;
output.ymin = ymin;  
output.bb_eval = bb_eval;
output.iteration = iteration;
output.nb_success = nb_success;
    
output.psize = psize;
output.psize_success = psize_success;
output.psize_max = psize_max;

if options.store_cache
    output.cache = cache;
else
    output.cache = "The option ''store_cache'' was set to false";
end

