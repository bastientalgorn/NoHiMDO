function Impossible_display(PB,x,v,w)

q = get_q(PB,x);
q2 = q./PB.q_scale;
qmin = -v./(2*w.^2);
q2min = qmin./PB.q_scale;

z1 = x(1:2);
z2 = x(3:4);
% Targets
t1 = q2min+z2;
t2 = z1-q2min;

lb = min(PB.lb);
ub = max(PB.ub);
K = PB.UserData.K;

hold off;
plot([-lb +lb],[0 0],'k');
hold on;
plot([0 0],[-lb +lb],'k');

xplot = [lb ub];
yplot = xplot+K;
plot(xplot,yplot,'--k');

xplot = 1 + exp(2*1i*pi*linspace(0,1,161));
plot(xplot,'--k');


plot(x(1),x(2),'or');
plot(x(3),x(4),'*b');
plot([z1(1) t1(1)],[z1(2) t1(2)],'--r');
plot([z2(1) t2(1)],[z2(2) t2(2)],'--b');



%axis([0.1 0.5 0.5 0.9])
axis equal

drawnow

