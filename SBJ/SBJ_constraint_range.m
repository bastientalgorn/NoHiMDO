function G4=SBJ_constraint_range(h,M,We,Wf,LD,SFC,Ws)

%-----THIS SECTION COMPUTES THE A/C RANGE (Breguet)-----%

if h<36089
     theta=1-0.000006875*Zh;
else
     theta=.7519;
end

Wt = We+Wf+Ws;

%-----THIS SECTION COMPUTES THE A/C RANGE-----%
range=(M*LD*661*sqrt(theta)/SFC)*log(Wt/(Wt-Wf));


% Range linearisation
% tau = M*661*sqrt(theta);
% LDref = 7.03;
% SFCref = 1;
% Wtref = 34300;
% Wfref = 10300;
% cLD = tau/SFCref*log(Wtref/(Wtref-Wfref));
% cSFC = -tau*LDref/(SFCref^2)*log(Wtref/(Wtref-Wfref));
% cWt = tau*LDref*SFCref*(1/Wtref - 1/(Wtref-Wfref));
% cWf = tau*LDref*SFCref/(Wtref-Wfref);
% R0 = (tau*LDref/SFCref)*log(Wtref/(Wtref-Wfref));
% range = R0+cLD*LD+cSFC*SFC+cWt*Wt+cWf*Wf;


G4= -range/2000+1;
% if G4>=0
%     disp(num2str([G4 Wt Wt-Wf LD/SFC],3));
% end
%global c_score 
%c_score = [c_score ; [G4 Wt Wt-Wf LD/SFC]];
