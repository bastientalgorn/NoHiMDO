%-------------------------------------------------------------------------------------%
%  NoHiMDO                                                                            %
%                                                                             f        %
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

function PB = build_problem(PB)
%close all
%clear all
%PB = Basic_problem_definition;

%===========================================
% If already built, then leave.
%===========================================
if isfield(PB,'is_built') && PB.is_built
    return;
end



%===========================================
% Read information from var
%===========================================
disp('================ Building problem =============');

% Number of variables (some of them can be vectors)
NV = length(PB.var);

% Name of the variable
PB.var_names = cell(0);
% Index of the subproblem in which this variable is DV or CV
PB.SP_index = [];
% Boolean that indicates if the variable is a CV
PB.is_coupling_variable    = [];
% Indexes of the variables to which the variable is linked
PB.links  = cell(0);
% Dimension of the variable
PB.var_dim = [];

% Number of components of X
NX = 0;
% Lower and upper bound for each components of x
PB.lb  = [];
PB.ub  = [];

% x_index_to_var_index: For each component x[i] of x, indicate what is the corresponding variable
PB.x_index_to_var_index = [];
% x_index_to_var_subindex: indicate what is the index of x[i] in this variable.
PB.x_index_to_var_subindex = [];

% Create a struct to have a fast "name to index" access.
PB.name_to_index_struct = struct;

% Loop on the variables
for i=1:NV
    v = PB.var{i};

    % Variable name
    var_name = v{1};
    if ~ischar(var_name)
        disp(['variable name is of type ' class(sp_index) ':']);
        disp(var_name)
        error(['Variable #' num2str(i) ': 1st argument (variable name) must be a char.']);
    end
    if strcmp(var_name,'auto')
        error('A variable cannot be named "auto". This is a reserved word');
    end
    PB.var_names{i} = var_name;
    PB.name_to_index_struct.(var_name) = i;
    

    % Subproblem index
    % Indicate in which subproblem the variable is involved
    sp_index = v{2};
    if ~isnumeric(sp_index) || numel(sp_index)~=1 || round(sp_index)~=sp_index
        disp(['subproblem index is of type ' class(sp_index) ':']);
        disp(sp_index);
        error(['Variable ' PB.var_names{i} ': 2nd argument (subproblem index) must be a scalar integer.']);
    end
    PB.SP_index(i) = int8(sp_index);

    % Coupling variable flag
    % Indicate if the variable is a design or coupling variable.
    % TRUE if the variable is a coupling variable
    is_cv = v{3};
    if (~isnumeric(is_cv) && ~islogical(is_cv)) || numel(is_cv)~=1 || logical(is_cv)~=is_cv
        disp(['coupling variable flag is of type ' class(is_cv) ':']);
        disp(is_cv);
        error(['Variable ' PB.var_names{i} ': 3rd argument (coupling variable flag) must be a scalar boolean.']);
    end
    PB.is_coupling_variable(i) = logical(is_cv);

    % Links (will be checked later...)
    PB.links{i}    = v{4};

    % dimension of the current variable
    dim = v{5};
    if ~isnumeric(dim) || numel(dim)~=1 || round(dim)~=dim
        disp(['variable dimension is of type ' class(dim) ':']);
        disp(dim);
        error(['Variable ' PB.var_names{i} ': 5th argument (variable dimension) must be a scalar integer.']);
    end
    PB.var_dim(i) = dim;

    % Lower bound
    lbi = v{6};
    % Check type and dimension
    if ~isnumeric(lbi) || min(size(lbi))>1
        disp(['lower bound is of type ' class(lbi) ':']);
        disp(lbi);
        error(['Variable ' PB.var_names{i} ': 6th argument (lower bound) must be a vector or scalar.']);
    end
    lb = zeros(1,dim);
    if length(lbi)==1
        lb(:) = lbi;
    elseif length(lbi)==dim
        lb = lbi(:)';
    else
        error(['Variable ' PB.var_names{i} ': 6th argument (lower bound) must be of dimension 1 or of the same dimension as the variable']);
    end
    PB.lb = [PB.lb , lb];

    % Upper bound
    ubi = v{7};
    % Check type and dimension
    if ~isnumeric(ubi) || min(size(ubi))>1
        disp(['upper bound is of type ' class(ubi) ':']);
        disp(ubi);
        error(['Variable ' PB.var_names{i} ': 7th argument (upper bound) must be a vector or scalar.']);
    end
    ub = zeros(1,dim);
    if length(ubi)==1
        ub(:) = ubi;
    elseif length(ubi)==dim
        ub = ubi(:)';
    else
        error(['Variable ' PB.var_names{i} ': 7th argument (upper bound) must be of dimension 1 or of the same dimension as the variable']);
    end
    PB.ub = [PB.ub , ub];

    % var_index_to_x_indexes : for a variable i, gives the corresponding indexes in x
    PB.var_index_to_x_indexes{i} = NX+(1:dim);
    % Total number of components of X
    NX = NX + dim;

    % Map allowing to know what each component of X represents.
    PB.x_index_to_var_index = [ PB.x_index_to_var_index , i*ones(1,dim) ];
    PB.x_index_to_var_subindex = [ PB.x_index_to_var_subindex , (1:dim) ];
end



% ================================================ %
% Analysis function
% ================================================ %

% Remove ".m" from analysis file
file = PB.analysis_file;
if length(file)>2 && strcmp(file(end-1:end),'.m')
    file = file(end-2);
end
% Check if analysis file exists
if exist(file,'file')~=2
    error(['Could not find analysis file ' file]);
end
PB.analysis_handle = str2func(file);
disp(['Analysis function: ' file]);

% ================================================ %
% end_of_iter function
% ================================================ %

% Remove ".m" from analysis file
if isfield(PB,'end_of_iter_file')
    file = PB.end_of_iter_file;
    if length(file)>2 && strcmp(file(end-1:end),'.m')
        file = file(end-2);
    end
    % Check if analysis file exists
    if exist(file,'file')~=2
        error(['Could not find endf_of_iter file ' file]);
    end
    disp(['End-of-iter function: ' file]);
end

% ======================================================= %
% Number of components, variables, and subproblems
% ======================================================= %

% Number of variables (some of them can be vectors)
PB.NV = NV;
disp(['Number of variables: ' num2str(NV)]);

% Number of subproblems
NSP = max(PB.SP_index);
PB.NSP = NSP;
disp(['Number of subproblems: ' num2str(NSP)]);

% Number of components in X
PB.NX = NX;
disp(['Number of components in X: ' num2str(NX)]);


% =============================== %
% CHECK
% =============================== %

% Check that all the variable names are different 
for i1=1:NV
    name_i1 = PB.var_names{i1};
    for i2=i1+1:NV
        name_i2 = PB.var_names{i2};
        if strcmp(name_i1,name_i2)
            disp(['Name of variable ' num2str(i1) ' : ' name_i1]);
            disp(['Name of variable ' num2str(i2) ' : ' name_i2]);
            error('Several variables have the same name');
        end
    end
    PB.var_names{i1} = cleanSpaces(name_i1);
end

% ================================================ %
% LINKS
% ================================================ %
link_matrix = false(NX,NX);

% interpret links and store them (in boolean form) in the link_matrix
for i1 = 1:NV
    L = PB.links{i1};
    i1_x_indexes = PB.var_index_to_x_indexes{i1};

    % For convenience, we convert L into a cell array.
    if ischar(L)
        if strcmp(L,'auto')
            % If the links are just defined with "auto", then try to automatically find the 
            % matching variables.
            L = get_automatic_links(i1,PB);
        else
            % If L is a char, convert into a cell array.
            L = {L};
        end    
    elseif isnumeric(L)
        % If L is a numeric, convert into a cell array.
        Lcopy = cell(0,0);
        for j=1:length(L)
            Lcopy{j} = L(j);
        end
        L = Lcopy;
        clear Lcopy;
    end
    
    % Then loop on the components of L
    for j=1:length(L)
        Lj = L{j};
        % if Lj is char convert to x indexes
        if ischar(Lj)
            % Search for index subset.
            % Exemple, variable u is linked to v(3:5)
            parenthesis_ind = min(find(Lj=='(' | Lj=='['));
            if ~isempty(parenthesis_ind)
                Lj_name = Lj(1:parenthesis_ind-1);
                Lj_ind  = eval(Lj(parenthesis_ind:end));
                i2 = name2varindex(Lj_name,PB);
                % TODO: Check dimension on Lj_ind
                % The following line will lead to an error if dimension
                % does not match. 
                % TODO: explicit error message.
                i2_x_indexes = PB.var_index_to_x_indexes{i2}(Lj_ind);
                disp('Subset index: ');
                disp(['Var name: ' Lj_name]);
                disp(['Indexes: ' num2str(Lj_ind)]);
                disp(['Corresponding indexes in x: ' num2str(i2_x_indexes)]);
            else
                i2 = name2varindex(Lj,PB);
                i2_x_indexes = PB.var_index_to_x_indexes{i2};
            end   
        elseif isnumeric(Lj)
            i2_x_indexes = PB.var_index_to_x_indexes{Lj};
        else
            disp(Lj)
            % TODO: more detailled message.
            error('Unrecognized link. Must be an integer or a string.');
        end
                
        % Check dimensions of i1_x_indexes and i2_x_indexes
        length_1 = length(i1_x_indexes);
        length_2 = length(i2_x_indexes);
        if ( length_1==length_2 )
            for k = 1:length_1
                link_matrix(i1_x_indexes(k),i2_x_indexes(k)) = true;
            end
        elseif (length_1==1)
            for k = 1:length_2
                link_matrix(i1_x_indexes,i2_x_indexes(k)) = true;
            end
        elseif (length_2==1)
            for k = 1:length_1
                link_matrix(i1_x_indexes(k),i2_x_indexes) = true;
            end    
        else
            disp(Lj)
            % TODO: more detailled message.
            error('Link size do not match');
        end
          
    end
end
clear i1_x_indexes i2_x_indexes k L Lj
PB.link_matrix = link_matrix;


% Create for each subproblem the indexes of the components of X that are :
% - a design variable
% - a coupling variable
% - a variable
for sp=1:NSP
    % Indexes of design variables for each subproblem
    D = [];
    for i=find(PB.SP_index==sp & ~PB.is_coupling_variable)
        PB.XDV_indexes{i} = length(D)+(1:PB.var_dim(i));
        D = [D ; PB.var_index_to_x_indexes{i}(:)];
    end
    PB.D_indexes{sp} = D;

    % Indexes of coupling variables
    C = [];
    for i=find(PB.SP_index==sp & PB.is_coupling_variable)
        PB.XCV_indexes{i} = length(C)+(1:PB.var_dim(i));
        C = [C ; PB.var_index_to_x_indexes{i}(:)];
    end
    PB.C_indexes{sp} = C;

    % Indexes of variables
    V = [];
    for i=find(PB.SP_index==sp)
        V = [V ; PB.var_index_to_x_indexes{i}(:)];
    end
    PB.V_indexes{sp} = V;
end
clear D C V





% Check the value of the subproblem associated with each variable 
% and detect if some variables are "dummy"
err = false;
PB.is_var_dummy = false(1,NV);
PB.is_x_dummy = false(1,NX);
for i=1:NV
    j_i = PB.SP_index(i);
    if (round(j_i)~=j_i)
        err = true;
        disp(['Variable ' PB.var_names{i} ...
              ' (#' num2str(i) ') belongs to subsystem ' ...
              num2str(j_i)]);
        disp('Not an integer!');
    end
    if (j_i>NSP)
        err = true;
        disp(['Variable ' PB.var_names{i} ...
              ' (#' num2str(i) ') belongs to subsystem ' ...
              num2str(j_i)]);
        disp(['Bigger than the larger subsystem index! (' num2str(NSP) ')']);
    end
    if (j_i<=0)
        disp(['Variable ' PB.var_names{i} ...
              ' (#' num2str(i) ...
              ') is detected as a "dummy" variable.']);
        PB.is_var_dummy(i) = true;
        PB.is_x_dummy(PB.var_index_to_x_indexes{i}) = true;
    end
end
if err
    error('Definition of j_i is not valid');
end

%===================================%
% Propagate and clean links         %
%===================================%

% Convert link matrix to double
link_matrix = double(link_matrix)+eye(size(link_matrix));
% Propagate
link_matrix_old = 0*link_matrix;
while any(any(link_matrix_old - link_matrix))
    link_matrix_old = link_matrix;
    link_matrix = link_matrix + link_matrix' + link_matrix*link_matrix;
    link_matrix = double(link_matrix>0);
end
clear link_matrix_old;
% Convert back to logical
link_matrix = (link_matrix>0);

% CHECK that no dummy variable is linked to another dummy variable
err = false;
list_of_dummy = find(PB.is_x_dummy);
for i1=list_of_dummy
    for i2=list_of_dummy
        if link_matrix(i1,i2) && (i1~=i2)
            disp(['Link between dummy variables ' num2str(i1) ' and ' num2str(i2) ' are not allowed!']);
            err = true;
        end
    end
end
if err
    error('Some dummy variables are linked');
end

%=============================================================
% If a comp-x is linked to a dummy comp-x, then they don't need to be connected 
% with any non-dummy comp-x.
%=============================================================

for i=1:NX
    % if a non-dummy variable i is connected with a dummy
    if ~PB.is_x_dummy(i) && any(link_matrix(i,PB.is_x_dummy))
        % Then i must not be connected with any other non-dummy
        link_matrix(i,~PB.is_x_dummy) = false;
    end
end

%======================================
% Parallel Coordination Matrix 
%======================================
% Matrix xavg_matrix is used to build the target value of each component of x-indexes
% prior to parallel optimization of the subsystems.
xavg_matrix = link_matrix;
% Clean dummy
for i=list_of_dummy
    xavg_matrix(i,:) = 0;
    xavg_matrix(i,i) = 1 ;
end
% Normalize lines
xavg_matrix = bsxfun( @ldivide , sum(xavg_matrix,2) , xavg_matrix );
% Check that lines are normalized correctly
for i=1:NX
    if (abs(sum(xavg_matrix(i,:))-1)>1e-12)
        error('xavg_matrix is not correctly normalized.');
    end
end
PB.xavg_matrix = xavg_matrix;



% Remove diagonal terms from link matrix
link_matrix(1:NX+1:NX*NX) = false;

%=============================================================
% Compute the list of the variables that are linked
%=============================================================

% L is the list of all the pairs of variables that are linked together
% Q_indexes{j} is the list of indexes (in L) that are considered within
% sub-problem j.
L = [];
for j=1:NSP
    PB.Q_indexes{j} = [];
end
for i1=1:NX
    for i2=i1+1:NX
        if link_matrix(i1,i2)
            % Add this pair to L
            L = [L  [i1;i2] ];
            NQ = size(L,2);
            % Add this index to the relevent Q_indexes
            
            % Variable corresponding to these x-indexes
            v1 = PB.x_index_to_var_index(i1);
            v2 = PB.x_index_to_var_index(i2);
            
            % Corresponding subproblems
            j1 = PB.SP_index(v1);
            j2 = PB.SP_index(v2);
            
            % Note: the Q_indexes allow to take only into accounts
            % the components of q relevent to the current subproblem.
            % That's nice if one of the other components of q is very 
            % large and make the penalty function very large (which 
            % may generate numeric problems)
            % Actually, the Q_indexes are generally not necessary.
            if j1>0
                PB.Q_indexes{j1}(end+1) = NQ;
            end
            if j2>0
                PB.Q_indexes{j2}(end+1) = NQ;
            end

        end
    end
end
clear v1 v2 j1 j2 i1 i2
PB.L = L;
PB.NQ = size(L,2);

for j=1:NSP
    PB.Q_indexes{j} = unique(PB.Q_indexes{j});
end


% Scale (inverse of the diff between bounds for each links couple)
S = PB.ub-PB.lb;
S = 1./S;
S(isinf(S))=1;
S(isnan(S))=1;
S(S<=0)=1;
S = 0.5*(S(L(1,:))+S(L(2,:)));
PB.q_scale = S;


for i=1:NX
    v  = PB.x_index_to_var_index(i);
    vi = PB.x_index_to_var_subindex(i);
    x_names = PB.var_names{v};
    if PB.var_dim(v)>1
        x_names = [x_names '(' num2str(vi) ')'];
    end
    PB.x_names{i} = x_names;
end


% End of built
PB.is_built = true;

PB
