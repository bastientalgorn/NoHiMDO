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
tab_problem = {'biquad','geo','sac','sbj'};
%tab_problem = {'geo'};

tab_pbfullname = {
    'Bi-Quadratic',
    'Geometric prog.',
    'Simplified MDO',
    'Supersonic business jet'
    };

tab_NI = [  8  16  32 64 128 256 512 ];
tab_NO = [512 256 128 64  32  16   8 ];
tab_gamma = [0.0 0.2 0.4 0.6 0.8 1.0];
tab_beta = [1.4 1.8 2.2 2.6];

tab_ylim = {[10^-20 10^0],[10^-3 10^+3];  % Quad
            [10^-20 10^0],[10^-2 10^10];  % Geo
            [10^-20 10^2],[10^-2 10^0];  % SAC / Simplified MDO
            [10^-2  10^1],[10^-2 10^1]}; % SBJ
        
tab_y0   = [0.5  100;
            0.08 12e+8;
            50 0.4061;
            7   3.8];

avg_metric = 'median';

% DISPLAY option
output_dir = '../results/';
figure_position = [362 63 1056 712];
tex_options = {'fontsize',30,'interpreter','latex','fontname','times'};
gca_fontsize = 20;
markersize = 8;

% Measures
inc_label = '$\epsilon_q$';
obj_label = '$\epsilon_f$';
tab_measure = {'inc','obj'};
tab_measure_label = { inc_label , obj_label };

beta_label = '$\beta$';
gamma_label = '$\gamma$';
iter_label = 'Cumulative number of iterations';
ni_label = 'Number of inner loop iterations';
no_label = 'Number of outer loop iterations';
newline = char(10);

tab_algo = {'mads','madsx','fmincon'};





tab_algo_fancy = tab_algo;
for i=1:3
    switch tab_algo{i}
        case 'fmincon'
            tab_algo_fancy{i} = 'fmincon';
        case 'mads'
            tab_algo_fancy{i} = 'Mads';
        case 'madsx'
            tab_algo_fancy{i} = 'Mads$\Delta$';
        otherwise
            error('algo not recognized');
            
    end
end
% Fix gamma, beta

% num2str precision
nsp = 2;


optim_method = 'mads';


do = [3 4];
% 1: courbes
% 2: boxplot
% 3: bar3
% 4: latex tabs
% 5: pareto (Do not use)
% 6: algo compare (Do not use)
% 7: boxplot algo compare
% 8: explication bounce

if ismember(1,do)
    gamma = 0.4;
    beta = 2.2;


    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)

        pb_name = tab_problem{ip};

        % LOOP ON THE ALGO
        for ia = 1:1%length(tab_algo)
            algo = tab_algo{ia};

            % LOOP ON THE MEASURES
            for im=1:2
                measure = tab_measure{im};
                measure_label  =  tab_measure_label{im};

                figure('position',figure_position,'color','w');
                hold on;
                leg_txt = cell(0);
                leg_ptr = zeros(0);

                % LOOP ON THE NI
                for in = 1:length(tab_NI)
                    NI = tab_NI(in);
                    NO = tab_NO(in);
                    color = get_color(in,length(tab_NI));
                    marker = get_marker(in,length(tab_NI));
                    x = (0:NO)*NI;
                    y = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_' avg_metric '_curve']);
                    y = [tab_y0(ip,im) y'];
                    y = max(y,1e-20);
                    leg_ptr(end+1) = plot(x,y,[marker '-'],'color',color);
                    s = '';
                    if NO<100
                        s = '$\;\,$';
                    end
                    if NO<10
                        s = '$\;\;\,\,$';
                    end
                    leg_txt{in} = ['NO=' num2str(NO) ',' s ' NI=' num2str(NI)];
                end
                set(gca,'yscale','log');

                xlabel(iter_label,tex_options{:});
                ylabel(measure_label,tex_options{:});
                xlim([0 NI*NO]);
                ylim(tab_ylim{ip,im});
                set(gca,'xtick',(0:1/8:1)*NI*NO,'fontsize',gca_fontsize);
                if ip==1 && im==1
                    legend(leg_ptr,leg_txt,tex_options{:},'fontsize',24,'location','best');
                    set(legend,'fontsize',24);
                end
                figure_name = ['Curve_' measure '_' algo '_'  pb_name];
                export_fig([output_dir figure_name '.pdf'],'-pdf');
                close(gcf);

            end % end measure
        end % end algo

    end % problem
end % member





%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
%     BOX PLOT         %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,do)
    xaxisfontsize = 22;
    xaxisleft = 0.2;
    xaxis0 = 0.55;

    gca_fontsize = 22;
    tex_options = {'fontsize',50,'interpreter','latex','fontname','times'};

    gamma = 0.4;
    beta = 2.2;
    
        
    tab_ylim_boxplot = [ 0 2 ; 0 15 ; 0 1 ; 1 0 ];
    tab_ylim_boxplot = 10.^tab_ylim_boxplot;
    color = [1 1 1];


    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};

        % LOOP ON THE MEASURES
        for im=1:2
            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im};

            % LOOP ON THE ALGO
            for ia = 1%:2
                algo = tab_algo{ia};

                figure('position',figure_position,'color','w');
                subplot(1,1,1);
                set(gca,'yscale','log');
                hold on;

                % fetch data
                xticklabel = cell(size(tab_NI));
                yl = [+inf -inf];

                % LOOP ON THE NI & NO
                for in = 1:length(tab_NI)
                    NI = tab_NI(in);
                    NO = tab_NO(in);
                    xticklabel{in} = ['NI=' num2str(NI) ', NO=' num2str(NO)];
                    vmean= get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_mean']);
                    vmin = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_min']);
                    vmax = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_max']);
                    v25  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.25']);
                    v50  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.50']);
                    v75  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.75']);

                    vmean = max(vmean,1e-20);
                    vmin  = max(vmin ,1e-20);
                    vmax  = max(vmax ,1e-20);
                    v25   = max(v25  ,1e-20);
                    v50   = max(v50  ,1e-20);
                    v75   = max(v75  ,1e-20);

                    yl(1) = min(yl(1),vmin);
                    yl(2) = max(yl(2),vmax);
                    
                    myboxplot(in,vmean,vmin,vmax,v25,v50,v75,color)
                end

                set(gca,'yscale','log','fontsize',gca_fontsize,'yminorgrid','off');
                ylabel(measure_label,tex_options{:});

                % yl = get(gca,'ylim');
                % yl(1) = yl(1)/((yl(2)/yl(1))^0.1);
                % ylim(yl);
                
                % Round the upper bound
                %yl(2) = 10^ceil(log10(yl(2)));
                yl(2) = tab_ylim_boxplot(ip,im);
                % get the lower bound down a little
                yl(1) = yl(1)/((yl(2)/yl(1))^0.15);
                ylim(yl);

                xtick = 1:length(tab_NI);
                %set(gca,'xtick',(1:length(tab_NI)),'xticklabel',[]);
                set(gca,'xtick',[],'xticklabel',[]);

                dx = 0;
                m1 = (yl(2)/yl(1))^0.080;
                m2 = (yl(2)/yl(1))^0.035;              
                for in = 1:length(tab_NI)
                    str_lineNI = num2str(tab_NI(in));
                    str_lineNO = num2str(tab_NO(in));
                    if in ==1
                        str_lineNI = ['NI = ' str_lineNI];
                        str_lineNO = ['NO = ' str_lineNO];
                        dx = -0.125;
                    else
                        dx=0;
                    end
                    pp(1) = text(xtick(in)+dx,yl(1)*m1,str_lineNI);
                    pp(2) = text(xtick(in)+dx,yl(1)*m2,str_lineNO);
                    for ipp=1:2
                        set(pp(ipp),'fontsize',xaxisfontsize,'horizontalalignment','center','interpreter','latex');
                    end
                end
                
                xlim([xaxisleft length(tab_NI)+0.5])

                figure_name = ['Boxplot_' measure '_' algo '_' pb_name ];
                export_fig([output_dir figure_name '.pdf'],'-pdf');
                close(gcf);

            end % end loop algo
        end % end loop measure
    end % end loop problem
end % end ismember







%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
%         BAR 3        %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(3,do)

    NI = 64;
    NO = 64;

    FLIP = false;
    ytick = (1:length(tab_beta));
    yticklabel = cell(1,length(tab_beta));
    for ib = 1:length(tab_beta)
        beta = tab_beta(ib);
        yticklabel{ib}=num2str(beta);
    end
    xtick = (1:length(tab_gamma));
    xticklabel = cell(1,length(tab_gamma));
    for ig = 1:length(tab_gamma)
        gamma = tab_gamma(ig);
        xticklabel{ig}=num2str(gamma);
    end

    if FLIP
        yticklabel = yticklabel(end:-1:1);
    end

    m = zeros(length(tab_beta),length(tab_gamma));

    % BAR3
    gca_fontsize = 25;

    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};

        % LOOP ON THE MEASURES
        for im=1:2
            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im};

            % LOOP ON THE ALGO
            for ia = 1%:2
                algo = tab_algo{ia};

                % FETCH DATA
                % LOOP ON beta & gamma
                for ib = 1:length(tab_beta)
                    beta = tab_beta(ib);
                    for ig = 1:length(tab_gamma)
                        gamma = tab_gamma(ig);
                        m(ib,ig) = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_' avg_metric]);
                    end
                end

                if FLIP
                    m = flipud(m);
                end


                figure('position',figure_position,'color','w');
                subplot(1,1,1);
                hold on;

                my2Dbarplot(m);
                set(gca,'xtick',xtick,'xticklabel',xticklabel,'ytick',ytick,'yticklabel',yticklabel);
                ylabel(beta_label,tex_options{:});
                xlabel(gamma_label,tex_options{:});
                zlabel(measure_label,tex_options{:});
                set(gca,'view',[-134 52],'fontsize',gca_fontsize);
                xlim([0.5 6.5])
                ylim([0.5 4.5])

                figure_name = ['Bar3_' measure '_' algo '_' pb_name ];
                export_fig([output_dir figure_name '.pdf'],'-pdf');

                close all
            end % end loop algo


        end % end loop measure

    end % end loop problem

end % end ismemnber






%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
%      LATEX TAB       %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(4,do)
optim_method = 'mads';
    gamma = 0.4;
    beta = 2.2;


    
    % Plot the INC and OBJ depending on NI
    % LATEX TAB
    out_file = [output_dir 'Tab_ByNI.tex'];
    fid = fopen(out_file,'w');

    fwrite(fid,'\begin{tabular}{| l |');
    for i=1:length(tab_NI)
        fwrite(fid,' | c c ');
    end
    fwrite(fid,' |}');
    fwrite(fid,newline);

    fwrite(fid,'\hline');
    fwrite(fid,newline);

    fwrite(fid,'NI ');
    for i=1:length(tab_NI)
        %fwrite(fid,[' & ' num2str(tab_NI(i)) ' & XXX ']);
        fwrite(fid,[' & \multicolumn{2}{|c|}{' num2str(tab_NI(i)) '} ']);
    end
    fwrite(fid,' \\');
    fwrite(fid,newline);

    fwrite(fid,'NO ');
    for i=1:length(tab_NI)
        %    fwrite(fid,[' & ' num2str(tab_NO(i)) ' & XXX ']);
        fwrite(fid,[' & \multicolumn{2}{|c|}{' num2str(tab_NO(i)) '} ']);
    end
    fwrite(fid,' \\');
    fwrite(fid,newline);

    fwrite(fid,'\hline\hline');
    fwrite(fid,newline);


    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};
        fwrite(fid,pb_name);
        fwrite(fid,newline);

        % fetch data
        for in = 1:length(tab_NI)
            NI = tab_NI(in);
            NO = tab_NO(in);
            obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric]);
            inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric]);
            fwrite(fid,[' & ' myformat(inc,true) ' & ' myformat(obj,true) ' % NI,NO = ' num2str(NI) ',' num2str(NO) newline]);
        end
        fwrite(fid,'\\\hline');
        fwrite(fid,newline);
    end

    fwrite(fid,'\end{tabular}');
    fclose(fid);


    
    
    % Plot the INC and OBJ depending on NI
    % LATEX TAB
    out_file = [output_dir 'Tab_ByPb1.tex'];
    fid = fopen(out_file,'w');

    fwrite(fid,'\begin{tabular}{| c c |');
    for i=1:length(tab_problem)
        fwrite(fid,' | c c ');
    end
    fwrite(fid,' |}');
    fwrite(fid,newline);

    fwrite(fid,'\hline');
    fwrite(fid,newline);

    fwrite(fid,'NI & NO ');
    for i=1:length(tab_problem)
        fwrite(fid,[' & \multicolumn{2}{|c|}{' tab_pbfullname{i} '} ']);
    end
    fwrite(fid,[' \\' newline]);

    fwrite(fid,'  &  ');
    for i=1:length(tab_problem)
        fwrite(fid,[' & ' inc_label ' & ' obj_label ' ']);
    end
    fwrite(fid,[' \\' newline]);
    fwrite(fid,'\hline\hline');
    fwrite(fid,newline);


    % fetch data
    for in = 1:length(tab_NI)

        NI = tab_NI(in);
        NO = tab_NO(in);
        fwrite(fid,['%=======================================' newline]);
        fwrite(fid,[num2str(NI) ' & ' num2str(NO) '  % NI, NO ' newline]);

        for ip = 1:length(tab_problem)
            pb_name = tab_problem{ip};
            obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric]);
            inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric]);
            fwrite(fid,[' & ' myformat(inc,true) ' & ' myformat(obj,true) ' % ' pb_name newline]);
        end
        fwrite(fid,'\\\hline');
        fwrite(fid,newline);
    end
    fwrite(fid,['%=======================================' newline]);
    fwrite(fid,'\end{tabular}');
    fclose(fid);








    % Plot the INC and OBJ depending on NI
    % LATEX TAB
    out_file = [output_dir 'Tab_ByPb2.tex'];
    fid = fopen(out_file,'w');

    fwrite(fid,'\begin{tabular}{| c || c');
    for i=1:length(tab_problem)
        fwrite(fid,' | c');
    end
    fwrite(fid,' |}');
    fwrite(fid,newline);

    fwrite(fid,'\cline{3-6}');
    fwrite(fid,newline);

    % Problem names
    fwrite(fid,'\multicolumn{2}{c|}{}');
    fwrite(fid,newline);
    fwrite(fid,' & Bi-Quadratic  ');
    fwrite(fid,newline);
    fwrite(fid,' & Geometric  ');
    fwrite(fid,newline);
    fwrite(fid,' & Simplified');
    fwrite(fid,newline);
    fwrite(fid,' & Supersonic \\');
    fwrite(fid,newline);
    fwrite(fid,'\multicolumn{2}{c|}{}   ');
    fwrite(fid,newline);
    fwrite(fid,' &');
    fwrite(fid,newline);
    fwrite(fid,' & Prog.');
    fwrite(fid,newline);
    fwrite(fid,' & MDO');
    fwrite(fid,newline);
    fwrite(fid,' & Business Jet \\');
    fwrite(fid,newline);

    fwrite(fid,'\hline');
    fwrite(fid,newline);

    % fetch data
    for in = 1:length(tab_NI)

        NI = tab_NI(in);
        NO = tab_NO(in);
        fwrite(fid,['%=======================================' newline]);
        fwrite(fid,['\multirow{2}{*}{$NI=' num2str(NI) '$, $NO=' num2str(NO) '$} % NI, NO ' newline]);


        fwrite(fid,['&' inc_label newline]);
        for ip = 1:length(tab_problem)
            pb_name = tab_problem{ip};
            inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric]);
            fwrite(fid,['   & ' myformat(inc,true) ' % ' pb_name newline]);
        end
        fwrite(fid,'\\\cline{2-6}');
        fwrite(fid,newline);

        %    fwrite(fid,['& ' newline']);
        fwrite(fid,['&' obj_label newline]);
        for ip = 1:length(tab_problem)
            pb_name = tab_problem{ip};
            obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric]);
            fwrite(fid,['   & ' myformat(obj,true) ' % ' pb_name newline]);
        end
        fwrite(fid,'\\\hline');
        fwrite(fid,newline);
    end
    fwrite(fid,['%=======================================' newline]);
    fwrite(fid,'\end{tabular}');
    fclose(fid);



    NO = 64;
    NI = 64;


    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};
        out_file = [output_dir 'Tab_ByBetaGamma_' pb_name '.tex'];
        fid = fopen(out_file,'w');
        fwrite(fid,'\begin{tabular}{| c | c || c ');
        for i=1:length(tab_gamma)
            fwrite(fid,' | c');
        end
        fwrite(fid,' |}');
        fwrite(fid,newline);
        fwrite(fid,['\cline{4-' num2str(3+length(tab_gamma)) '}' newline]);

        fwrite(fid,['\multicolumn{3}{c|}{}' newline]);
        fwrite(fid,['& \multicolumn{' num2str(length(tab_gamma)) '}{|c|}{ $\gamma$ }\\']);
        fwrite(fid,newline);
        fwrite(fid,['\cline{4-' num2str(3+length(tab_gamma)) '}' newline]);

        fwrite(fid,['\multicolumn{3}{c|}{}' newline]);
        for g=tab_gamma
            fwrite(fid,['& ' num2str(g) ' ']);
        end
        fwrite(fid,[' \\' newline]);
        fwrite(fid,['\hline' newline]);

        fwrite(fid,[newline '% -------------------------------' newline newline]);

        nr = num2str(2*length(tab_beta));
        fwrite(fid,['\multirow{' nr '}{*}{$\beta$}' newline]);

        for ib=1:length(tab_beta)
            beta = tab_beta(ib);
            %fwrite(fid,['& ' num2str(beta)]);
            fwrite(fid,[newline '% -------------------------------' newline newline]);
            fwrite(fid,['% beta = ' num2str(beta) newline]);
            fwrite(fid,['& \multirow{2}{*}{' num2str(beta) '}' newline]);
            fwrite(fid,['    % inc:' newline]);
            fwrite(fid,['    & ' inc_label ' ' newline '   ']);
            for ig=1:length(tab_gamma)
                gamma = tab_gamma(ig);
                inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric]);
                fwrite(fid,[' & ' myformat(inc,true)]);
            end
            fwrite(fid,['\\' newline]);

            fwrite(fid,['& %empty' newline]);

            fwrite(fid,['    % obj:' newline]);
            fwrite(fid,['    & ' obj_label ' ' newline '   ']);
            for ig=1:length(tab_gamma)
                gamma = tab_gamma(ig);
                obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric]);
                fwrite(fid,[' & ' myformat(obj,true)]);
            end
            fwrite(fid,['\\' newline]);
            if ib<length(tab_beta)
                fwrite(fid,['\cline{2-' num2str(3+length(tab_gamma)) '}' newline]);
            else
                fwrite(fid,['\hline' newline]);
            end
        end

        fwrite(fid,[newline '% -------------------------------' newline newline]);
        fwrite(fid,'\end{tabular}');
        fclose(fid);


    end






    NO = 64;
    NI = 64;

    out_file = [output_dir 'Tab_ByBetaGamma_AllPb.tex'];
    fid = fopen(out_file,'w');
    fwrite(fid,'\begin{tabular}{| c | c || c ');
    for i=1:length(tab_gamma)
        fwrite(fid,' | c');
    end
    fwrite(fid,' |}');
    fwrite(fid,newline);

    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};
        pb_fullname = tab_pbfullname{ip};

        if ip>1
            fwrite(fid,['\multicolumn{' num2str(3+length(tab_gamma)) '}{c}{}\\[0pt]']);
            fwrite(fid,newline);
        end

        fwrite(fid,['\multicolumn{3}{c}{}' newline]);
        fwrite(fid,['& \multicolumn{' num2str(length(tab_gamma)) '}{c}{ {\normalsize ' pb_fullname '} }\\']);
        fwrite(fid,newline);


        fwrite(fid,['\cline{4-' num2str(3+length(tab_gamma)) '}' newline]);

        fwrite(fid,['\multicolumn{3}{c|}{}' newline]);
        fwrite(fid,['& \multicolumn{' num2str(length(tab_gamma)) '}{|c|}{ $\gamma$ }\\']);
        fwrite(fid,newline);
        fwrite(fid,['\cline{4-' num2str(3+length(tab_gamma)) '}' newline]);

        fwrite(fid,['\multicolumn{3}{c|}{}' newline]);
        for g=tab_gamma
            fwrite(fid,['& ' num2str(g) ' ']);
        end
        fwrite(fid,[' \\' newline]);
        fwrite(fid,['\hline' newline]);

        fwrite(fid,[newline '% -------------------------------' newline newline]);

        nr = num2str(2*length(tab_beta));
        fwrite(fid,['\multirow{' nr '}{*}{$\beta$}' newline]);

        for ib=1:length(tab_beta)
            beta = tab_beta(ib);
            %fwrite(fid,['& ' num2str(beta)]);
            fwrite(fid,[newline '% -------------------------------' newline newline]);
            fwrite(fid,['% beta = ' num2str(beta) newline]);
            fwrite(fid,['& \multirow{2}{*}{' num2str(beta) '}' newline]);
            fwrite(fid,['    % inc:' newline]);
            fwrite(fid,['    & ' inc_label ' ' newline '   ']);
            for ig=1:length(tab_gamma)
                gamma = tab_gamma(ig);
                inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric]);
                fwrite(fid,[' & ' myformat(inc,true)]);
            end
            fwrite(fid,['\\' newline]);

            fwrite(fid,['& %empty' newline]);

            fwrite(fid,['    % obj:' newline]);
            fwrite(fid,['    & ' obj_label ' ' newline '   ']);
            for ig=1:length(tab_gamma)
                gamma = tab_gamma(ig);
                obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric]);
                fwrite(fid,[' & ' myformat(obj,true)]);
            end
            fwrite(fid,['\\' newline]);
            if ib<length(tab_beta)
                fwrite(fid,['\cline{2-' num2str(3+length(tab_gamma)) '}' newline]);
            else
                fwrite(fid,['\hline' newline]);
            end
        end

        fwrite(fid,[newline '% -------------------------------' newline newline]);

    end

    fwrite(fid,'\end{tabular}');
    fclose(fid);

end






% ========== %
%   PARETO   %
% ========== %
if ismember(5,do)

    gamma = 0.4;
    beta = 2.2;

    % Plot the INC curves
    for ip = 1:1%length(tab_problem)
        pb_name = tab_problem{ip};
        figure('position',figure_position,'color','w');
        hold on;
        leg_txt = cell(0);
        leg_ptr = zeros(0);

        for in = 1:length(tab_NI)
            NI = tab_NI(in);
            NO = tab_NO(in);
            color = get_color(in,length(tab_NI));
            marker = get_marker(in,length(tab_NI));
            curve_inc = get_value(pb_name,NO,NI,gamma,beta,optim_method,['inc_' avg_metric '_curve']);
            curve_obj = get_value(pb_name,NO,NI,gamma,beta,optim_method,['obj_' avg_metric '_curve']);
            leg_ptr(in) = plot(curve_obj,curve_inc,[marker '-'],'color',color);
            s = '';
            if NO<100
                s = '$\;\,$';
            end
            if NO<10
                s = '$\;\;\,\,$';
            end
            leg_txt{in} = ['NO=' num2str(NO) ',' s ' NI=' num2str(NI)];
        end
        set(gca,'xscale','log','yscale','log');

        xlabel(obj_label,tex_options{:});
        ylabel(inc_label,tex_options{:});
        if ip==1
            legend(leg_ptr,leg_txt,tex_options{:},'location','best');
            set(legend,'fontsize',30);
        end
        %xlim([0 4096]);
        %set(gca,'xtick',(0:512:4096),'fontsize',gca_fontsize);
        figure_name = ['IncVsObj_' pb_name];
        export_fig([output_dir figure_name '.pdf'],'-pdf');
        close(gcf);

    end
end



% =================== %
%   ALGO COMPARISON   %
% =================== %
if ismember(6,do)

    nb_optim_method = 3;
    gamma = 0.4;
    beta = 2.2;
    tab_NI = [  8  16  32 64 128 256 512 ];
    tab_NO = [512 256 128 64  32  16   8 ];

    % LOOP ON PROBLEMS
    for ip = 1:length(tab_problem)
        disp('=====================')
        pb_name = tab_problem{ip}

        
        % LOOP ON MEASURES
        for im = 1:2

            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im}
                
            figure('position',figure_position,'color','w');
            hold on;
            leg_txt = cell(0);
            leg_ptr = zeros(0,0);

            
            % LOOP ON NI & NO
            for in = 1:3%length(tab_NI)
                NI = tab_NI(in);
                NO = tab_NO(in);


                for io = 1:nb_optim_method
                    switch io
                        case 1
                            optim_method = 'fmincon';
                            leg_txt{io} = 'fmincon';
                            marker = '-';
                        case 2
                            optim_method = 'mads';
                            leg_txt{io} = 'MADS';
                            marker = 'o-';
                        case 3
                            optim_method = 'madsx';
                            leg_txt{io} = 'MADS-$\Delta$';
                            marker = '*-';
                        otherwise
                            io
                            nb_optim_method
                            error('"in" is not a proper optim method index');
                    end
                    color = get_color(in,length(tab_NI));
                    x = (0:NO)*NI;
                    y = get_value(pb_name,NO,NI,gamma,beta,optim_method,[measure '_q0.50_curve']);
                    y = [tab_y0(ip,im) y'];
                    %leg_ptr(io) = plotwithmarkers(x,y,10,'marker',marker,'color',color,'linestyle','-','nicemarkers','markersize',markersize,'offset',io/nb_optim_method);
                    leg_ptr(io) = plot(x,y,marker,'color',color);
                    % 25th and 75th percentiles
                    %y = get_value(pb_name,NO,NI,gamma,beta,optim_method,[measure '_q0.25_curve']);
                    %plot(x,y,[marker '--'],'color',color);
                    %y = get_value(pb_name,NO,NI,gamma,beta,optim_method,[measure '_q0.75_curve'])
                    %plot(x,y,[marker '--'],'color',color);
                end
                set(gca,'yscale','log');

            end
            
            xlabel(iter_label,tex_options{:});
            ylabel(measure_label,tex_options{:});
            if ip==1 && im==1
                legend(leg_ptr,leg_txt,tex_options{:},'location','best');
                set(legend,'fontsize',30);
            end
            xlim([0 NI*NO]);
            ylim(tab_ylim{ip,im});
            set(gca,'xtick',(0:1/8:1)*NI*NO,'fontsize',gca_fontsize);
            figure_name = ['AlgoCompare_' measure '_' pb_name ];
            export_fig([output_dir figure_name '.pdf'],'-pdf');
            %pause
            close(gcf);

        end

    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%     Algo comparison with BoxPlot         %
%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(7,do)
    xaxisfontsize = 22;
    xaxisleft = -0.2;
    xaxis0 = 1.25;

    gca_fontsize = 22;
    tex_options = {'fontsize',50,'interpreter','latex','fontname','times'};

    gamma = 0.4;
    beta = 2.2;
    
    color = { [1 1 1],0.85*[1 1 1],0.5*[1 1 1] };
    
    tab_ylim_boxplot = [ 0 2 ; 0 15 ; 5 1 ; 1 1 ];
    tab_ylim_boxplot = 10.^tab_ylim_boxplot;
    
    % LOOP ON THE PROBLEMS
    for ip = 1:length(tab_problem)
        pb_name = tab_problem{ip};

        % LOOP ON THE MEASURES
        for im=1:2
            measure = tab_measure{im};
            measure_label  =  tab_measure_label{im};

            figure('position',figure_position,'color','w');
            subplot(1,1,1);
            set(gca,'yscale','log');
            hold on;
            
            k = 1;
            yl = [+inf -inf];
            ytext = -inf;
            % LOOP ON THE NI & NO
            for in = 1:length(tab_NI)
                NI = tab_NI(in);
                NO = tab_NO(in);
                %xticklabel{in} = ['NI=' num2str(NI) ', NO=' num2str(NO)];
                
                % LOOP ON THE ALGO
                for ia = 1:3
                    algo = tab_algo{ia};

                    % fetch data
                    xticklabel = cell(size(tab_NI));

                    vmean= get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_mean']);
                    vmin = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_min']);
                    vmax = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_max']);
                    v25  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.25']);
                    v50  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.50']);
                    v75  = get_value(pb_name,NO,NI,gamma,beta,algo,[measure '_q0.75']);
                    vmean = max(vmean,1e-20);
                    vmin  = max(vmin ,1e-20);
                    vmax  = max(vmax ,1e-20);
                    v25   = max(v25  ,1e-20);
                    v50   = max(v50  ,1e-20);
                    v75   = max(v75  ,1e-20);
                    
                    myboxplot(k,vmean,vmin,vmax,v25,v50,v75,color{ia})
                    
                   
                    yl(1) = min(yl(1),vmin);
                    yl(2) = max(yl(2),vmax);
                    k = k+1;
                    if in==1
                        ytext = max(ytext,vmax);
                    end
            
                end

                k = k+1;
            end % end loop algo
            
    
            set(gca,'yscale','log','fontsize',gca_fontsize);
            ylabel(measure_label,tex_options{:});

            % Round the upper bound
            %yl(2) = 10^ceil(log10(yl(2)));
            yl(2) = tab_ylim_boxplot(ip,im);
            % get the lower bound down a little
            yl(1) = yl(1)/((yl(2)/yl(1))^0.15);
            ylim(yl);

            xtick = 4*(1:length(tab_NI))-2;
            set(gca,'xtick',[],'xticklabel',[],'yminorgrid','off');

            dx = 0;
            m1 = (yl(2)/yl(1))^0.080;
            m2 = (yl(2)/yl(1))^0.035;
            for in = 1:length(tab_NI)
                str_lineNI = num2str(tab_NI(in));
                str_lineNO = num2str(tab_NO(in));
                if in ==1
                    str_lineNI = ['NI = ' str_lineNI];
                    str_lineNO = ['NO = ' str_lineNO];
                    dx = -0.5;
                else
                    dx=0;
                end
                pp(1) = text(xtick(in)+dx,yl(1)*m1,str_lineNI);
                pp(2) = text(xtick(in)+dx,yl(1)*m2,str_lineNO);
                for ipp=1:2
                    set(pp(ipp),'fontsize',xaxisfontsize,'horizontalalignment','center','interpreter','latex');
                end
            end

            for in = 1:length(tab_NI)
                plot(4*in*[1 1],ylim,'--','color',0.5*[1 1 1]);
            end
            xlim([-1 4*length(tab_NI)+0.01])
            
            if ip==1
                for ia=1:3
                    pp = text(ia-0.15,ytext*2,tab_algo_fancy{ia});
                    set(pp,'fontsize',xaxisfontsize,...
                           'horizontalalignment','left',...
                           'verticalalignment','middle',...
                           'interpreter','latex','rotation',90);
                end
            end
            
            figure_name = ['Boxplot_Compare_' measure '_' pb_name ];
            export_fig([output_dir figure_name '.pdf'],'-pdf');
            close(gcf);


        end % end loop measure
    end % end loop problem
end % end ismember



%========%
% BOUNCE %
%========%
if ismember(8,do)
    gamma = 0.4;
    beta = 2.2;
    in  = 4;
    NI = tab_NI(in);
    NO = tab_NO(in);
    color = get_color(in,length(tab_NI));
    marker = get_marker(in,length(tab_NI));
    ip =1;
    pb_name = tab_problem{ip};
    algo = 'mads';
    
    x = (0:NO)*NI;
    eq = get_value(pb_name,NO,NI,gamma,beta,algo,['inc_median_curve']);
    eq = [tab_y0(ip,1) eq'];
    ef = get_value(pb_name,NO,NI,gamma,beta,algo,['obj_median_curve']);
    ef = [tab_y0(ip,2) ef'];


    figure;
    subplot(3,1,1);
    plot(x,eq,[marker '-'],'color',color);
    xlabel(iter_label,tex_options{:});
    ylabel(tab_measure_label{1},tex_options{:});
    set(gca,'yscale','log');
    xlim([0 512]);
    set(gca,'xtick',(0:1/8:1)*NI*NO,'fontsize',gca_fontsize);

    subplot(3,1,2);
    plot(x,ef,[marker '-'],'color',color);
    xlabel(iter_label,tex_options{:});
    ylabel(tab_measure_label{2},tex_options{:});
    xlim([0 512]);
    set(gca,'xtick',(0:1/8:1)*NI*NO,'fontsize',gca_fontsize);



    subplot(3,1,3);
    f = ef+3;
    plot(x,f,[marker '-'],'color',color);
    xlabel(iter_label,tex_options{:});
    ylabel(tab_measure_label{2},tex_options{:});
    xlim([0 512]);
    set(gca,'xtick',(0:1/8:1)*NI*NO,'fontsize',gca_fontsize);


    figure_name = ['Bounce_'  pb_name];
    export_fig([output_dir figure_name '.pdf'],'-pdf');
    %close(gcf);


end % member



