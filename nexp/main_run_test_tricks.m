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
close all
clear all

% ================== %
%    PARAM LISTS     %
% ================== %

global NEXP

TEST = 0;


switch TEST
    case 0
        tab_pb = {'geo','sac','sbj','gsr'};
        tab_optionString = {'','1','2','3','4','5','6','7','8'};        
        NEXP = 40;
        NI = 100;
        NO = 100;
    case 1
        tab_pb = {'geo'};
        %tab_optionString = {'','1','2','3'};
        tab_optionString = {'5'};
        NEXP = 8;
        NI = 5;
        NO = 10;
    case 2
        tab_pb = {'sbj'};
        tab_optionString = {'','1'};
        NEXP = 2;
        NI = 50;
        NO = 50;
        
end





default_run_options.display = false;
default_run_options.oracle = false;
default_run_options.filter = false;
default_run_options.reduced_bounds = false;
default_run_options.w_scheme = 'normal';
default_run_options.cache = false;
default_run_options.constraints_cv = false;
default_run_options.realistic_obj = false;
% Algo stops if this inconsistency is reached
default_run_options.inc_stop = 1e-6;
% Algo stops if the inconsistency has not decreased for that many iter
default_run_options.noprogress_stop = 20;
% Mads-Delta : special initialization of the poll size
default_run_options.mads_delta = false;
% Number of Inner/Outer loop iterations
default_run_options.NI = NI;
default_run_options.NO = NO;
% Hyper-parameters of the penalty update scheme
default_run_options.beta = 2.0;
default_run_options.gamma = 0.5;
% Initial value for the w vector
default_run_options.w0 = 1;








tab_options = cell(size(tab_optionString));
for i=1:length(tab_optionString)
    tab_optionString{i} = unique(sort(tab_optionString{i}));
    use_k = ismember('123456789',tab_optionString{i});
    run_options = default_run_options;
    for k=1:9
        use = use_k(k);
        if (use_k(k))
            switch k
                case 1
                    run_options.cache = use;
                case 2
                    run_options.constraints_cv = use;
                case 3
                    run_options.realistic_obj = use;
                case 4
                    run_options.reduced_bounds = use;
                case 5
                    run_options.filter = use;
                case 6
                    run_options.oracle = use;
                case 7
                    if use
                        run_options.w_scheme = 'median';
                    end
                case 8
                    if use
                        run_options.w_scheme = 'max';
                        if use_k(7)
                            error('Option 7 and 8 cannot be used at the same time');
                        end
                    end
                case 9
                    error('Not defined');
                otherwise
                    error('Not defined');
            end

        end
    end
    tab_options{i} = run_options;
end

% Build list of runs
list_run = [];
for i=1:length(tab_pb)
    for j=1:length(tab_optionString)
        for k=1:NEXP
            list_run(end+1,:) = [i j k];
        end
    end
end
N = length(list_run);
list_run = list_run(randperm(N),:);

disp('===========================================================');
for r=1:N
    i = list_run(r,1);
    j = list_run(r,2);
    k = list_run(r,3);
    pb = tab_pb{i};
    optionString = ['no' num2str(NO) '_ni' num2str(NI) '_s' tab_optionString{j}];
    options = tab_options{j};
    nexp = k;
    try
        run_nexp(pb,options,optionString,nexp);
    catch exception
        display_exception(exception);
    end
end
disp('===========================================================');
%exit;


