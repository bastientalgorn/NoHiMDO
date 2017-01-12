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

function run_nexp(pb_name,options,optionString,nexp)

dir = './matfiles/';
! mkdir -p ./matfiles
% =================================
% LOAD PROBLEM
% =================================
PB = define_problem(pb_name);
%fs = get_f_star(pb_name);
x0 = load_x0(PB);
global NEXP
%NEXP = size(x0,1);
%NEXP = 5;
if (nexp>NEXP)
    error('nexp > NEXP');
end

% =================================
% RUN
% =================================
out_file = [dir pb_name '_' optionString ];

disp(['out_file : ' out_file]);
if exist([out_file '.mat'],'file')
    disp('EXISTS! (skip)');
    return;
else

    disp(['run_nexp: building ' out_file]);
    out_nexp_file = [out_file '_nexp' num2str(nexp) '.mat'];
    disp(['    out_nexp_file : ' num2str(nexp)]);
    % Check if the mat file for this nexp is existing
    if exist(out_nexp_file,'file')
        load(out_nexp_file);
        disp('    Exists! (skip)');
    else
        % If not, do the math
        try
            [inc,obj] = NoHiSolver(PB,x0(nexp,:),options);
            save(out_nexp_file,'inc','obj');
        catch exception
            display_exception(exception);
        end
        % Check if the file exists. If not, create a crash file
        % (Just so that the user knows that this given calculation
        % crashed.
        if ~exist(out_nexp_file,'file')
            disp('Create crash file');
            system(['touch ' out_nexp_file '.crash']);
        else
            % If the file exist (ie, the calculation worked fine)
            % , we remove the eventual crash file (which might have
            % been created during a previous run).
            system(['rm ' out_nexp_file '.crash 2>/dev/null']);
        end

    end

end

build = true;
for nexp=1:NEXP
    out_nexp_file = [out_file '_nexp' num2str(nexp) '.mat'];
    if ~exist(out_nexp_file,'file')
        build = false;
        return;
    end
end


if build
    tab_inc = [];
    tab_obj = [];
    for nexp=1:NEXP
        out_nexp_file = [out_file '_nexp' num2str(nexp) '.mat'];
        load(out_nexp_file);
        tab_inc = [tab_inc , inc(:)];
        tab_obj = [tab_obj , obj(:)];
    end
    
    obj_doublemed = zeros(size(tab_inc,1),1);
    for i=1:size(tab_inc,1)
        inc_i = tab_inc(i,:);
        inc_med = median(inc_i);
        obj_i = tab_obj(i,:);
        obj_i = obj_i(inc_i<=inc_med);
        obj_doublemed(i) = median(obj_i);
    end

    obj_q3 = zeros(size(tab_inc,1),1);
    for i=1:size(tab_inc,1)
        inc_i = tab_inc(i,:)<1e-3;
        obj_i = tab_obj(i,:);
        obj_i(~inc_i) = 1;
        obj_q3(i) = median(obj_i);
    end
    obj_q6 = zeros(size(tab_inc,1),1);
    for i=1:size(tab_inc,1)
        inc_i = tab_inc(i,:)<1e-6;
        obj_i = tab_obj(i,:);
        obj_i(~inc_i) = 1;
        obj_q6(i) = median(obj_i);
    end


    mat = [ mean(tab_inc') ,
            std(tab_inc') ,
            min(tab_inc') ,
            max(tab_inc') ,
            quantile(tab_inc',0.10) ,
            quantile(tab_inc',0.25) ,
            quantile(tab_inc',0.50) ,
            quantile(tab_inc',0.75) ,
            quantile(tab_inc',0.90) ,
            mean(tab_obj') ,
            std(tab_obj') ,
            min(tab_obj') ,
            max(tab_obj') ,
            quantile(tab_obj',0.10) ,
            quantile(tab_obj',0.25) ,
            quantile(tab_obj',0.50) ,
            quantile(tab_obj',0.75) ,
            quantile(tab_obj',0.90) ,
            obj_doublemed',
            obj_q3',
            obj_q6' ];
    mat = mat';
    save([out_file '.mat'],'mat');
    disp(out_file)
end


