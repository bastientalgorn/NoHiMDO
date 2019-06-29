function [obj,y,c_ineq] = SBJ_subsystem_analysis(subproblem_index,x_DV,PB)

% default values
obj = 0;
y = [];
c_ineq = [];

% Constants
h = PB.UserData.altitude;
M = PB.UserData.Mach;

switch subproblem_index
    case 1
        [SFC,We,LD,Ws,Wf] = get_variable(x_DV,PB,'SFC_1','We_1','LD_1','Ws_1','Wf_1');  
        Wt = SBJ_obj_range(We,Wf,Ws);
        y = Wt;
        obj = Wt;
        c_ineq = SBJ_constraint_range(h,M,We,Wf,LD,SFC,Ws);
        c_ineq = max(c_ineq,0)^2;
    case 2
        [D,T] = get_variable(x_DV,PB,'D_2','T');  
        ineqtol = get_variable(x_DV,PB,'ineqtol');  
        [SFC,We,ESF] = SBJ_obj_power(h,M,T,D);
        y = [We,SFC,ESF];
        obj = 0;
        c_ineq = SBJ_constraint_power(h,M,T,D);%-ineqtol;
    case 3
        [ESF,Wt,theta,tc,ARw,LAMBDAw,Sref,Sht,ARht,LAMBDAht,Lw,Lht] = split_vector(x_DV);
        Z = [tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht];
        [L,D,LD] = SBJ_obj_dragpolar(Z,LAMBDAht,Lw,Lht,Wt,theta,ESF);
        y = [D,LD,L];
        obj = 0;
        c_ineq = SBJ_constraint_dragpolar(Z,Lw,Lht,Wt,tc);
    case 4
        [L,tc,ARw,LAMBDAw,Sref,Sht,ARht,lambda] = get_variable(x_DV,PB,'L_4','tc_4','ARw_4','LAMBDAw_4','Sref_4','Sht_4','ARht_4','lambda_4');  
        [t,ts] = get_variable(x_DV,PB,'t','ts');  

        Z = [tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht];
        [Ws,Wf,theta] = SBJ_obj_weight(Z,t,ts,lambda,L);
        y = [Ws,Wf,theta];
        obj = 0;
        c_ineq = SBJ_constraint_weight(Z,t,ts,lambda,L);
    otherwise
        error('unrecognized subproblem index')
end