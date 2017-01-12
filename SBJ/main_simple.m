clear all
clc


%----INITIALIZE ALL VARIABLES, SYSTEM AND LOCAL, AND UPPER/LOWER BOUNDS----(problem dependent)%


[tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht] = split_vector([0.05  55000.0  1.6  3.0  60.0  500.0  100.0 5.5]);

%------t---------------ts-------------tpr_r-Sw_ht-Lw----Lht--T
%X_lb=[0.1*ones(1,9)   0.1*ones(1,9)  0.1   40.0  0.01  1.0  0.1];
X=    [3.0*ones(1,9)   6.0*ones(1,9)  0.3   45.0  0.15  1.5  0.6];

t  = X(1:9);
ts = X(10:18);
[lambda,LAMBDAht,Lw,Lht,T] = split_vector(X(19:23));
clear X
%X_ub=[4.0*ones(1,9)   9.0*ones(1,9)  0.4   70.0   0.2  3.5  1.0];

%                                   twst

[L,We,Wt,theta,ESF,D,Wf,LD,SFC,Ws] = split_vector([ 25000.0  15000.0 25000.0  10.0    1.0    40000.0 25000   5.0    2.0   25000]);


clear X Y Z

%==================================================================================



% SP 4
Z = [tc,nan,nan,ARw,LAMBDAw,Sref,Sht,ARht];
[Ws , Wf, theta] = SBJ_obj_weight(Z,t,ts,lambda,L);
Z = [tc,nan,nan,ARw,LAMBDAw,Sref,nan,nan];
G1               = SBJ_constraint_weight(Z,t,ts,lambda,L);

% SP 3
Z = [tc,h,M,ARw,LAMBDAw,Sref,Sht,ARht];
[L,D,LD] = SBJ_obj_dragpolar(Z,LAMBDAht,Lw,Lht,Wt,theta,ESF);
G2       = SBJ_constraint_dragpolar(Z,Lw,Lht,Wt,tc);


% SP 2
[SFC,We,ESF] = SBJ_obj_power(h,M,T,D);
G3           = SBJ_constraint_power(h,M,T,D);

% SP 1
Wt = SBJ_obj_range(We,Wf,Ws   );
G4 = SBJ_constraint_range(h,M,We,Wf,LD,SFC,Ws   );



%---------L(1)-----We(1)-----Wt(4)----twst(2)----ESF(2)----D(3)----Wf(4)----L/D(4)----SFC(4)----Ws(1)
Y_out=[   L        We        Wt       theta      ESF       D       Wf       LD        SFC       Ws];

%----Combine constraint vectors----%
G=[G1 G2 G3 G4];


Y = [L,We,Wt,theta,ESF,D,Wf,LD,SFC,Ws];
Y_in = Y;






disp('=============================================');
disp('=============================================');
disp('=============================================');
Y_in_2 = Y_in;
Y_out_2 = Y_out;
G_2 = G;

load check.mat
if ~isequal(Y_in,Y_in_2)
    error('NON CONSISTENT Y_in');
end
if ~isequal(Y_out,Y_out_2)
    error('NON CONSISTENT Y_out');
end
if ~isequal(G,G_2)
    error('NON CONSISTENT G');
end

dispOK;