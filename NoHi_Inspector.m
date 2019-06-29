% Inspector !!

close all
clear all
load output.mat
PB = build_problem(PB);
whos

addpath('./plotmethods');

output

NOmax = size(output.subproblems,2);
NSP = PB.NSP;
figure; hold on;
for nsp=1:NSP
    disp(nsp);
    D = PB.D_indexes{nsp};
    C = PB.C_indexes{nsp};
    
    for no=1:NOmax
        s = output.subproblems(nsp,no);
        x_Dm = s.x_master(D);
        x_Cm = s.x_master(C);
        x_Dt = s.x_target(D);
        x_Ct = s.x_target(C);
        x_Df = s.x_D;
        x_Cf = s.x_C;
        
        
        % Compute xf: value of x at the end of subproblem.
        xf = s.x_master;
        xf(D) = s.x_D;
        xf(C) = s.x_C;
        
        qm = get_q(PB,s.x_master);
        qf = get_q(PB,xf);
        
        qmaxm = max(abs(qm));
        qmaxf = max(abs(qf));
        
        
        color = get_color(nsp,NSP);
        subplot(2,2,1); hold on;
        plot(no+[0 nsp/(NSP+1)],[qmaxm qmaxf],'-','color',color);
        %plot(no+0.5*(2*nsp-1)/NSP,qmaxf,'v','color',get_color(nsp,NSP));
        plot(2:length(output.tab_inc)+1,output.tab_inc,'.k');
        set(gca,'yscale','log');
        plot(no,s.psize(nsp),'o','color',color);
        ylabel('q and psize');
        xlabel('no');
        
        subplot(2,2,2); hold on;
        plot(no,qmaxf/qmaxm,'o','color',color);
        xlabel('no');
        ylabel('qmaster/qfinal');
        
        subplot(2,2,3); hold on;
        plot(qmaxf/qmaxm,s.psize(nsp),'o','color',color);
        xlabel('qmaster/qfinal');
        ylabel('psize');
        set(gca,'xscale','log','yscale','log');
        
    end
end
