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

function [obj,y,c_ineq] = subsystem_analysis(PB,j,x_DV)

obj = 0;
y = [];
c_ineq = [];

switch PB.name
    
    case 'PbClass1'
        switch j
            case 1
                % ATC - Subproblem 1
                [fu,fa] = split_vector(x_DV);
                obj = fu+fa;
                y = [];
            case 2
                % ATC - Subproblem 2
                L = 0.1;
                [a,b,c] = split_vector(x_DV);
                fa = a.^3+1/(b+L)+L*sqrt(c);
                w = log(1+fa+c);
                y = [w,fa];
                obj = 0;
            case 3
                % ATC - Subproblem 3
                L = 0.1;
                [u,v,w] = split_vector(x_DV);
                fu = sqrt(u+v)+1/(v+L)+L*w^2;
                c = log(11-w)+fu+1/(u+L);
                y = [c,fu];
                obj = 0;
                
            otherwise
                error('unrecognized subproblem index')
        end
        
    case 'PbClass2'
        switch j
            case 1
                % ATC - Subproblem 1
                L = 0.1;
                [u,v,w] = split_vector(x_DV);
                fu = sqrt(u+v)+1/(v+L)+L*w^2;
                c = log(11-w)+fu+1/(u+L);
                y = [c,fu];
                c_ineq = v-3*u;
                obj = 0;
            case 2
                % ATC - Subproblem 2
                L = 0.1;
                [a,b,c] = split_vector(x_DV);
                fa = a.^3+1/(b+L)+L*sqrt(c);
                w = log(1+fa+c);
                y = [w,fa];
                obj = 0;
            case 3
                % ATC - Subproblem 3
                [fu,fa] = split_vector(x_DV);
                obj = fu+fa;
                y = [];
            otherwise
                error('unrecognized subproblem index')
        end
        
    case 'PbClass3'
         switch j
            case 1
                % ATC - Subproblem 1
                L = 0.1;
                [u,v,w] = split_vector(x_DV);
                fu = sqrt(u+v)+1/(v+L)+L*w^2;
                c = log(11-w)+fu+1/(u+L);
                y = [c,fu];
                c_ineq = v-3*u;
                obj = 0;
            case 2
                % ATC - Subproblem 2
                L = 0.1;
                [a,b,c] = split_vector(x_DV);
                fa = a.^3+1/(b+L)+L*sqrt(c);
                w = log(1+fa+c);
                y = [w,fa];
                obj = 0;
            case 3
                % ATC - Subproblem 3
                [fu,fa] = split_vector(x_DV);
                obj = fu+fa;
                h = fu-6*fa;
                y = [h];
            otherwise
                error('unrecognized subproblem index')
        end
        
        
    case 'biquad'
        switch j
            case -1
                % For Optimization All-At-Once
                obj = (x_DV-1)^2+(x_DV+1)^2;
            case 1
                % ATC - Subproblem 1
                obj = x_DV(1)+x_DV(2);
            case 2
                % ATC - Subproblem 2
                a = (x_DV-1)^2;
                obj = 0;
                y = a;
            case 3
                % ATC - Subproblem 3
                b = (x_DV+1)^2;
                obj = 0;
                y = b;
            otherwise
                error('unrecognized subproblem index')
        end

    case 'geo'
        switch j
            case 1
                % ATC - Subproblem 1
                [z3,z4,z5,z6,z7] = split_vector(x_DV);
                z1 = sqrt(z3^2+z4^-2+z5^2);
                z2 = sqrt(z5^2+z6^2 +z7^2);
                obj = z1^2+z2^2;
                c_ineq(1) = z3^-2 + z4^+2 - z5^+2;         
                c_ineq(2) = z5^+2 + z6^-2 - z7^+2;
            case 2
                % ATC - Subproblem 2
                [z8,z9,z10,z11] = split_vector(x_DV);
                z3 = sqrt(z8^2+z9^-2+z10^-2+z11^2);
                y = z3;
                obj = 0;
                c_ineq(1) = z8^+2 + z9^+2  - z11^+2;
                c_ineq(2) = z8^-2 + z10^+2 - z11^+2;
            case 3
                % ATC - Subproblem 3
                [z11,z12,z13,z14] = split_vector(x_DV);
                z6 = sqrt(z11^2+z12^2+z13^2+z14^2);
                y = z6;
                obj = 0;
                c_ineq(1) = z11^+2 + z12^-2 - z13^+2;
                c_ineq(2) = z11^+2 + z12^+2 - z14^+2;
            otherwise
                error('unrecognized subproblem index')
        end


    case 'gsr'
        switch j
            case 1
              [x1,x2,x3] = split_vector(x_DV);
              % SP objective
              F1 = 0.7854*x1*x2^2*(3.3333*x3*x3 + 14.9335*x3 - 43.0934);
              y = F1;
              % Constraints
              g5 = 27/(x1*x2^2*x3) -1;
              g6 = 397.5/(x1*x2^2*x3^2) -1;
              g9 = x2*x3/40 -1;
              g10 = 5*x2/x1 -1;
              g11 = x1/(12*x2) -1;
              c_ineq = [g5,g6,g9,g10,g11];
            case 2
              [x1,x2,x3,x4,x6] = split_vector(x_DV);
              % SP objective
              F2 = -1.5079*x1*x6^2;
              F4 = 7.477*x6^3;
              F6 = 0.7854*x4*x6^2;
              y = F2+F4+F6;
              % Constraints
              g1 = sqrt( ((745*x4)/(x2*x3))^2 + 1.69e+7)/(110*x6^3) -1;
              g3 = (1.5*x6 + 1.9)/x4 -1;
              g7 = 1.93*x4^3/(x2*x3*x6^4) -1;
              c_ineq = [g1,g3,g7];
            case 3
              [x1,x2,x3,x5,x7] = split_vector(x_DV);
              % SP objective
              F3 = -1.5079*x1*x7^2;
              F5 = 7.477*x7^3;
              F7 = 0.7854*x5*x7^2;
              y = F3+F5+F7;
              % Constraints
              g2 = sqrt( ((745*x5)/(x2*x3))^2 + 1.575e+8)/(85*x7^3) -1;
              g4 = (1.1*x7 + 1.9)/x5 -1;
              g8 = 1.93*x5^3/(x2*x3*x7^4) -1;
              c_ineq = [g2,g4,g8];
            case 4
              [f1,f2,f3] = split_vector(x_DV);
              obj = f1+f2+f3;

            otherwise
                error('unrecognized subproblem index')
        end



    case 'sac'
        switch j
            case 1
                % ATC - Subproblem 1
                [Ww,Wf] = split_vector(x_DV);
                obj = 60000 + Ww + Wf;
            case 2
                % ATC - Subproblem 2
                [x01,x02,x21,x22] = split_vector(x_DV);
                x0 = [x01,x02];
                x2 = [x21,x22];
                Ww = 4000*(1+norm(x0-1)^2)*(1+norm(x2-1)^2);
                y = Ww;
                obj = 0;
            case 3
                % ATC - Subproblem 3
                [x01,x02,x31,x32] = split_vector(x_DV);
                x0 = [x01,x02];
                x3 = [x31,x32];
                xs = 10*(x0+x3);
                eggholder_handle = @(x) -(x(2)+47) * sin(sqrt(abs(x(2)+x(1)/2+47))) - x(1) * sin(sqrt(abs(x(1)-(x(2)+47))));
                EH = eggholder_handle(xs);
                omega = (1+norm(x0-2)^2)*(1+0.001*norm(x3-2)^2)*(1+1000*abs(EH));
                Dr = 0.025+0.004*log10(omega);
                Wf = 20000 + 380952*Dr + 9523809*Dr*Dr;
                y = Wf;
                obj = 0;
            otherwise
                error('unrecognized subproblem index')
        end


    case 'sbj'
        h = 55000;
        M = 1.4;
        switch j
            case 1
                [SFC,We,LD,Ws,Wf] = split_vector(x_DV);
                Wt = SBJ_obj_range(We,Wf,Ws);
                y = Wt;
                obj = Wt;
                c_ineq = SBJ_constraint_range(h,M,We,Wf,LD,SFC,Ws);
                c_ineq = max(c_ineq,0)^2;
            case 2
                [D,T] = split_vector(x_DV);
                [SFC,We,ESF] = SBJ_obj_power(h,M,T,D);
                y = [We,SFC,ESF];
                obj = 0;
                c_ineq = SBJ_constraint_power(h,M,T,D);
            case 3
                [ESF,Wt,theta,tc,ARw,LAMBDAw,Sref,Sht,ARht,LAMBDAht,Lw,Lht] = split_vector(x_DV);
                Z = [tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht];
                [L,D,LD] = SBJ_obj_dragpolar(Z,LAMBDAht,Lw,Lht,Wt,theta,ESF);
                y = [D,LD,L];
                obj = 0;
                c_ineq = SBJ_constraint_dragpolar(Z,Lw,Lht,Wt,tc);
            case 4
                [L,tc,ARw,LAMBDAw,Sref,Sht,ARht,lambda] = split_vector(x_DV(1:8));
                t  = x_DV(9:17);
                ts = x_DV(18:26);
                Z = [tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht];
                [Ws,Wf,theta] = SBJ_obj_weight(Z,t,ts,lambda,L);
                y = [Ws,Wf,theta];
                obj = 0;
                c_ineq = SBJ_constraint_weight(Z,t,ts,lambda,L);
            otherwise
                error('unrecognized subproblem index')
        end

            
            
end


