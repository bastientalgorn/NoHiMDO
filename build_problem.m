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

function PB = build_problem(PB)

disp('================ Building problem =============');

% Number of variables (some of them can be vectors)
NV = length(PB.var);

% Name of the variable
PB.varnames = cell(0);
% Index of the subproblem in which this variable is DV or CV
PB.SP_index = [];
% Boolean that indicates if the variable is a CV
PB.cv    = [];
% Indexes of the variables to which the variable is linked
PB.links  = cell(0);
% Dimension of the variable
PB.dim = [];

% Number of components of X
NX = 0;
% Lower and upper bound for each components of x
PB.lb  = [];
PB.ub  = [];

% For each component x[i] of x, indicate (left col) what is the corresponding
% variable, and (right col) what is the index of x[i] in this variable.
PB.x2v = [];

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
    PB.varnames{i} = var_name;

    % Subproblem index
    % Indicate in which subproblem the variable is involved
    sp_index = v{2};
    if ~isnumeric(sp_index) || numel(sp_index)~=1 || round(sp_index)~=sp_index
        disp(['subproblem index is of type ' class(sp_index) ':']);
        disp(sp_index);
        error(['Variable ' PB.varnames{i} ': 2nd argument (subproblem index) must be a scalar integer.']);
    end
    PB.SP_index(i) = int8(sp_index);

    % Coupling variable flag
    % Indicate if the variable is a design or coupling variable.
    % TRUE if the variable is a coupling variable
    is_cv = v{3};
    if (~isnumeric(is_cv) && ~islogical(is_cv)) || numel(is_cv)~=1 || logical(is_cv)~=is_cv
        disp(['coupling variable flag is of type ' class(is_cv) ':']);
        disp(is_cv);
        error(['Variable ' PB.varnames{i} ': 3rd argument (coupling variable flag) must be a scalar boolean.']);
    end
    PB.cv(i) = logical(is_cv);

    % Links (will be checked later...)
    PB.links{i}    = v{4};

    % dimension of the current variable
    dim = v{5};
    if ~isnumeric(dim) || numel(dim)~=1 || round(dim)~=dim
        disp(['variable dimension is of type ' class(dim) ':']);
        disp(dim);
        error(['Variable ' PB.varnames{i} ': 5th argument (variable dimension) must be a scalar integer.']);
    end
    PB.dim(i) = dim;


    % Lower bound
    lbi = v{6};
    % Check type and dimension
    if ~isnumeric(lbi) || min(size(lbi))>1
        disp(['lower bound is of type ' class(lbi) ':']);
        disp(lbi);
        error(['Variable ' PB.varnames{i} ': 6th argument (lower bound) must be a vector or scalar.']);
    end
    lb = zeros(1,dim);
    if length(lbi)==1
        lb(:) = lbi;
    elseif length(lbi)==dim
        lb = lbi(:)';
    else
        error(['Variable ' PB.varnames{i} ': 6th argument (lower bound) must be of dimension 1 or of the same dimension as the variable']);
    end
    PB.lb = [PB.lb , lb];

    % Upper bound
    ubi = v{7};
    % Check type and dimension
    if ~isnumeric(ubi) || min(size(ubi))>1
        disp(['upper bound is of type ' class(ubi) ':']);
        disp(ubi);
        error(['Variable ' PB.varnames{i} ': 7th argument (upper bound) must be a vector or scalar.']);
    end
    ub = zeros(1,dim);
    if length(ubi)==1
        ub(:) = ubi;
    elseif length(ubi)==dim
        ub = ubi(:)';
    else
        error(['Variable ' PB.varnames{i} ': 7th argument (upper bound) must be of dimension 1 or of the same dimension as the variable']);
    end
    PB.ub = [PB.ub , ub];

    % X_indexes : indexes of the components of the variable in X
    PB.X_indexes{i} = NX+(1:dim);
    % Total number of components of X
    NX = NX + dim;

    % Map allowing to know what each component of X represents.
    PB.x2v = [ PB.x2v ; [ i*ones(dim,1) , (1:dim)' ] ];

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
disp(['Analysis function: ' file]);


% ======================================================= %
% Number of components, variables, subproblems and links
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

% Check that all the variable names are different and do not contain
% spaces
for i1=1:NV
    name_i1 = PB.varnames{i1};
    for i2=i1+1:NV
        name_i2 = PB.varnames{i2};
        if strcmp(name_i1,name_i2)
            disp(['Name of variable ' num2str(i1) ' : ' name_i1]);
            disp(['Name of variable ' num2str(i2) ' : ' name_i2]);
            error('Several variables have the same name');
        end
    end
    PB.varnames{i1} = cleanSpaces(name_i1);
end

% ================================================ %
% LINKS
% ================================================ %

% Replace strings or cells by corresponding integer.
for i = 1:NV
    L = PB.links{i};
    if strcmp(PB.varnames{i},L)
        error(['Variable ' num2str(i) ' (' L ') is linked to itself.']);
    end
    if ischar(L)
        % Find which variable name matches L
        i2 = name2varindex(L,PB)
        PB.links{i} = i2;
    elseif iscell(L)
        L_int = zeros(1,length(L));
        % Loop on the content of L
        for n=1:length(L)
            Ln = L{n};
            % Find which variable name matches Ln
            i2 = name2varindex(Ln,PB);
            L_int(n) = i2;
        end
        PB.links{i} = L_int;
    else
        % Store values as a row array
        PB.links{i} = L(:)';
    end
end

% Display links.
for i1=1:NV
    for i2 = PB.links{i1}
        disp([ PB.varnames{i1} ' (#' num2str(i1)  ') <------> '...
            PB.varnames{i2} ' (#' num2str(i2) ')']);
    end
end


% Create for each subproblem the indexes of the components of X that are :
% - a design variable
% - a coupling variable
% - a variable
for sp=1:NSP
    % Indexes of design variables for each subproblem
    D = [];
    for i=find(PB.SP_index==sp & ~PB.cv)
        PB.XDV_indexes{i} = length(D)+(1:PB.dim(i));
        D = [D ; PB.X_indexes{i}(:)];
    end
    PB.D_indexes{sp} = D;


    % Indexes of coupling variables
    C = [];
    for i=find(PB.SP_index==sp & PB.cv)
        PB.XCV_indexes{i} = length(C)+(1:PB.dim(i));
        C = [C ; PB.X_indexes{i}(:)];
    end
    PB.C_indexes{sp} = C;

    % Indexes of variables
    V = [];
    for i=find(PB.SP_index==sp)
        V = [V ; PB.X_indexes{i}(:)];
    end
    PB.V_indexes{sp} = V;
end
clear D C V


% Check that no variable is a target to itself
% And that the links are sym.
TR_check = false(NV,NV);
for i=1:NV
    L_i = PB.links{i};
    for link=L_i
        if i==link
            i
            link
            L_i
            error('i == link');
        end
        TR_check(i,link) = true;
    end
end
err = false;
for i1=1:NV-1
    for i2=i1+1:NV
        if TR_check(i1,i2)~=TR_check(i2,i1)
            disp(['Link of variables ' num2str(i1) ...
                ' and ' num2str(i2) ' are not symmetric']);
            err = true;
        end
    end
end
if err
    error('Some links are not symmetric');
end


% Check values of j_i and detect if some variables are "dummy"
err = false;
PB.dummy = false(1,NV);
for i=1:NV
    j_i = PB.SP_index(i);
    if (round(j_i)~=j_i)
        err = true;
        disp(['Variable ' PB.varnames{i} ' (#' num2str(i) ') belongs to subsystem ' num2str(j_i)]);
        disp('Not an integer!');
    end
    if (j_i>NSP)
        err = true;
        disp(['Variable ' PB.varnames{i} ' (#' num2str(i) ') belongs to subsystem ' num2str(j_i)]);
        disp(['Bigger than the larger subsystem index! (' num2str(NSP) ')']);
    end
    if (j_i<=0)
        disp(['Variable ' PB.varnames{i} ' (#' num2str(i) ') is a "dummy" variable.']);
        PB.dummy(i) = true;
    end
end
if err
    error('Definition of j_i is not valid');
end



% Compute the list of the variables that are linked
L = [];
for j=1:NSP
    PB.Q_indexes{j} = [];
end
for i1=1:NV
    j1 = PB.SP_index(i1);
    if ~isempty(PB.links{i1})
        for i2=PB.links{i1}
            j2 = PB.SP_index(i2);
            if i2>i1
                l1 = size(L,2)+1;
                L = [L  [PB.X_indexes{i1};PB.X_indexes{i2}] ];
                l2 = size(L,2);
                if j1>0
                    PB.Q_indexes{j1} = [ PB.Q_indexes{j1} ; (l1:l2)' ];
                end
                if j2>0
                    PB.Q_indexes{j2} = [ PB.Q_indexes{j2} ; (l1:l2)' ];
                end
            end
        end
    end
end
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


% Display variable names of coupled variables
k = 0;
for i=1:NV % Loop on the variables
    name_i = PB.varnames{i};
    if PB.cv(i)
        name_i = [name_i ' (CV #'];
    else
        name_i = [name_i ' (DV #'];
    end
    name_i = [name_i num2str(i) ', dim ' num2str(PB.dim(i)) ')'];

    for link=PB.links{i} % What are the TR for this variable?
        if i<link
            name_link = PB.varnames{link};
            if PB.cv(link)
                name_link = [name_link ' (CV #'];
            else
                name_link = [name_link ' (DV #'];
            end
            name_link = [name_link num2str(link) ', dim ' num2str(PB.dim(link)) ')'];
            s = ['Link ' num2str(k) ': ' name_i ' -> ' name_link ];
            disp(s);
            k = k+1;
        end
    end
end

% Check that all links are symmetrical
for i=1:NV % Loop on the variables
    for link=PB.links{i} % What are the TR for this variable?
        if ~ismember(i,PB.links{link})
            i
            link
            error('Link not symmetrical');
        end
    end
end
disp('Test link symmetry : OK');

% Check that linked variables have equal bounds
% and equal dimension
err_bounds = false;
err_dim    = false;
for i1=1:NV % Loop on the variables
    for i2=PB.links{i1} % What are the TR for this variable?

        % Check bounds for linked variables that are not dummy
        if i1<i2 && ~PB.dummy(i1) && ~PB.dummy(i2)
            X_indexes1 = PB.X_indexes{i1};
            X_indexes2 = PB.X_indexes{i2};

            lb1 = PB.lb(X_indexes1);
            lb2 = PB.lb(X_indexes2);
            if ~isequal(lb1,lb2)
                disp(['Variable ' num2str(i1) ' (' PB.varnames{i1} ') and ' num2str(i2)...
                    ' (' PB.varnames{i2} ') does not have equal lower bound.']);
                disp(['lb(' num2str(i1) ') = ' num2str(PB.lb(i1)) '; lb(' num2str(i2) ') = ' num2str(PB.lb(i2))]);
                err_bounds = false;
            end

            ub1 = PB.ub(X_indexes1);
            ub2 = PB.ub(X_indexes2);
            if ~isequal(ub1,ub2)
                disp(['Variable ' num2str(i1) ' (' PB.varnames{i1} ') and ' num2str(i2)...
                    ' (' PB.varnames{i2} ') does not have equal upper bound.']);
                disp(['ub(' num2str(i1) ') = ' num2str(PB.ub(i1)) '; ub(' num2str(i2) ') = ' num2str(PB.ub(i2))]);
                err_bounds = true;
            end
        end

        % Check dimension for all linked variables
        if i1<i2
            dim1 = PB.dim(i1);
            dim2 = PB.dim(i2);
            if dim1~=dim2
                disp(['Variable ' num2str(i1) ' (' PB.varnames{i1} ') and ' num2str(i2)...
                    ' (' PB.varnames{i2} ') does not have equal dimension.']);
                disp(['ub(' num2str(i1) ') = ' num2str(PB.dim(i)) '; ub(' num2str(i2) ') = ' num2str(PB.dim(i2))]);
                err_dim = true;
            end
        end
    end
end
if err_bounds
    error('Two linked variables should have the same bounds (except if one of them is dummy)');
end
if err_dim
    error('Two linked variables should have the same dimension');
end
disp('Test bound equality : OK');

% End of built
PB.is_built = true;

% Check
PB
x1 = rand(NX,1);
for i=1:NV
    if ~PB.cv(i)
        disp('====================');
        disp(PB.varnames{i});
        sp = PB.SP_index(i)
        XDV = x1(PB.D_indexes{sp})
        disp(['XDV indexes: ' num2str(PB.XDV_indexes{i})]);
        XDVi = XDV(PB.XDV_indexes{i})
        Vi = x1(PB.X_indexes{i})
        if ~isequal(XDVi,Vi)
            XDVi
            Vi
            error('Not equal...');
        end
    end
end



