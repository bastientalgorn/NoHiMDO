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
disp('============ PLOT RESULTS ============================');

addpath('./plotmethods');

skip = true;

% Parameters lists
% tab_problem = {'biquad','geo','sac','sbj'};
% tab_problem = {'biquad','geo'};
% tab_problem = {'geo','sac','gsr','sbj'};
% tab_algo = {'','1','2','3'};
% NI = 5;
% NO = 10;

tab_problem = {'sbj'};
tab_problem = {'geo','sac','gsr','sbj'};
tab_algo = {'','1','2','3','4','5','6','7'};
NI=100;
NO=100;

EPSILON = 1e-6;



tab_ylim = {[10^-6 10^-3],[10^-6 10^+6];  % Geo
            [10^-3 10^0],[10^-3 10^0];  % Sac
            [10^-3 10^0],[10^-3 10^0];  % gsf
            [10^-3 10^0],[10^-3 10^0]}; % sbj

tab_ylim = {[10^-6 10^-4],[10^-3 10^+5];  % Geo
            [10^-6 10^+0],[10^-2 10^+0];  % Sac
            [10^-6 10^+0],[10^-5 10^+0];  % gsf
            [10^-3 10^+0],[10^-4 10^+0]}; % sbj
        
% tab_y0 = [  1 0 ;
%     1 0;
%     1 0;
%     1 0];

tab_y0 = ones(4,2);

avg_metric = 'q3';

% DISPLAY option
output_dir = './results/';
figure_position = [362 63 1056 712];
tex_options = {'fontsize',30,'interpreter','latex','fontname','times'};
gca_fontsize = 20;
markersize = 8;

% Measures
inc_label = '$\epsilon_q$';
obj_label = '$\epsilon_f$';
tab_measure = {'inc','obj'};
tab_measure_label = { inc_label , obj_label };

iter_label = 'Cumulative number of iterations';
ni_label = 'Number of inner loop iterations';
no_label = 'Number of outer loop iterations';
newline = char(10);




x = (0:NO)*NI;
xtick = (0:NO)*NI;
if NO==100
    xtick = (0:20:NO)*NI;
end


tab_algo_fancy = tab_algo;
for i=1:length(tab_algo)
    switch tab_algo{i}
        case ''
            tab_algo_fancy{i} = 'Default';
        case '1'
            tab_algo_fancy{i} = 'Cache (1)';
        case '2'
            tab_algo_fancy{i} = 'Constraints CV (2)';
        case '3'
            tab_algo_fancy{i} = 'Realistic Obj. (3)';
        case '4'
            tab_algo_fancy{i} = 'Reduced bounds (4)';
        case '5'
            tab_algo_fancy{i} = 'Filter (5)';
        case '6'
            tab_algo_fancy{i} = 'Oracle (6)';
        case '7'
            tab_algo_fancy{i} = 'Median scheme (7)';
        otherwise
            error('algo not recognized');

    end
end

% num2str precision
nsp = 2;


do = [1];
% 1: courbes
% 2: boxplot
% 3: bar3
% 4: latex tabs
% 5: pareto (Do not use)
% 6: algo compare (Do not use)
% 7: boxplot algo compare
% 8: explication bounce

if ismember(1,do)

    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)

        pb_name = tab_problem{ip};

        % LOOP ON THE MEASURES
        for im=1:2
            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im};

            figure('position',figure_position,'color','w');
            hold on;
            leg_txt = cell(0);
            leg_ptr = zeros(0);


            % LOOP ON THE ALGO
            for ia = 1:length(tab_algo)
                algo = tab_algo{ia};


                color = get_color(ia,length(tab_algo));
                marker = get_marker(ia,length(tab_algo));

                y = get_value_trick(pb_name,NO,NI,algo,[measure '_' avg_metric '_curve']);
                y = [tab_y0(ip,im) y'];
                y = max(y,EPSILON);
                leg_ptr(end+1) = plot(x,y,[marker '-'],'color',color);
                leg_txt{ia} = tab_algo_fancy{ia};

            end % end algo

            set(gca,'yscale','log');

            xlabel(iter_label,tex_options{:});
            ylabel(measure_label,tex_options{:});
            xlim([0 NI*NO]);
            %ylim(tab_ylim{ip,im});
            set(gca,'xtick',xtick,'fontsize',gca_fontsize);
            if true %ip==1 && im==1
                legend(leg_ptr,leg_txt,tex_options{:},'fontsize',24,'location','best');
                set(legend,'fontsize',24);
            end
            figure_name = ['Curve_' measure '_'  pb_name];
            export_fig([output_dir figure_name '.pdf'],'-pdf');
            %kitkit
            %close(gcf);


        end % end measure


    end % problem
end % member





%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
%     BOX PLOT         %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,do)
    yaxisfontsize = 22;
    yaxisleft = 0.2;
    yaxis0 = 0.55;

    gca_fontsize = 22;
    tex_options = {'fontsize',50,'interpreter','latex','fontname','times'};

    gamma = 0.4;
    beta = 2.2;


    tab_xlim_boxplot = [ 0 2 ; 0 15 ; 0 1 ; 1 0 ];
    tab_xlim_boxplot = 10.^tab_xlim_boxplot;
    color = [1 1 1];


    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};
        figure('position',figure_position,'color','w');
        
        % LOOP ON THE MEASURES
        for im=1:2
            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im};

            subplot(1,2,im);
            if im==1
                set(gca,'xscale','log');
            end
            hold on;

 

            % LOOP ON THE ALGO
            for ia = 1:length(tab_algo)
                algo = tab_algo{ia};

                xl = [+inf -inf];

                vmean= get_value_trick(pb_name,NO,NI,algo,[measure '_mean']);
                vmin = get_value_trick(pb_name,NO,NI,algo,[measure '_min']);
                vmax = get_value_trick(pb_name,NO,NI,algo,[measure '_max']);
                v25  = get_value_trick(pb_name,NO,NI,algo,[measure '_q0.25']);
                v50  = get_value_trick(pb_name,NO,NI,algo,[measure '_q0.50']);
                v75  = get_value_trick(pb_name,NO,NI,algo,[measure '_q0.75']);

                vmean = max(vmean,EPSILON);
                vmin  = max(vmin ,EPSILON);
                vmax  = max(vmax ,EPSILON);
                v25   = max(v25  ,EPSILON);
                v50   = max(v50  ,EPSILON);
                v75   = max(v75  ,EPSILON);

                xl(1) = min(xl(1),vmin);
                xl(2) = max(xl(2),vmax);
                myboxplot(ia,vmean,vmin,vmax,v25,v50,v75,color,2)
            end

            set(gca,'xscale','log','fontsize',gca_fontsize,'xminorgrid','off');
            xlabel(measure_label,tex_options{:});

            

            if im==1
                v = xl(2)
                v = 10^ceil(ceil(v));
                xl(2) = v*1.5;
                xl(1) = EPSILON/2;
                xl;
                xlim(xl);
                r = log10(EPSILON):v;
                if length(r)>4
                    if mod(length(r),2)==1
                        r = r(1:2:end);
                    elseif mod(length(r),3)==1
                        r = r(1:3:end);
                    else
                        r = [r(1:2:end-1) r(end)+1];
                    end                
                end
                xtick = 10.^(r);
            else
                xl(2) = min(xl(2),1e+6);
            end
                
            set(gca,'xtick',xtick);

            ytick = 1:length(tab_algo);
            if im==1
                set(gca,'ytick',ytick,'yticklabel',tab_algo_fancy,'ydir','reverse');
            else
                set(gca,'ytick',ytick,'yticklabel',[],'ydir','reverse');
            end


            ylim([yaxisleft length(tab_algo)+0.5])

            figure_name = ['Boxplot_' pb_name ];
            export_fig([output_dir figure_name '.pdf'],'-pdf');


        end % end loop measure
        %close(gcf);
    end % end loop problem
end % end ismember






