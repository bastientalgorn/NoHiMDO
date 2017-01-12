function G1 = SBJ_constraint_weight(Z,t,ts,lambda,L)

%----Inputs----%
t=t(:)/12;%convert to feet
ts=ts(:)/12;%convert to feet

t1=t(1:3); 
t2=t(4:6);  
t3=t(7:9); 
ts1=ts(1:3); 
ts2=ts(4:6); 
ts3=ts(7:9); 
G=4000000*144;
E=10600000*144;
nu=0.3;

%----Inputs----%

beta=.9;
[c,c_box,Sweep_40,D_mx,b,a]=Wing_Mod(Z,lambda);
teq1=((t1.^3)/4+(3*t1).*(ts1-t1/2).^2).^(1/3);
teq2=((t2.^3)/4+(3*t2).*(ts2-t2/2).^2).^(1/3);
teq3=((t3.^3)/4+(3*t3).*(ts3-t3/2).^2).^(1/3);
l=0.6*c_box;
h=Z(1)*(beta*c(1:3))-(0.5)*(ts1+ts3);
A_top=(0.5)*(t1.*l)+(1/6)*(t2.*h);
A_bottom=(0.5)*(t3.*l)+(1/6)*(t2.*h);
Y_bar=h.*(2*A_top)./(2*A_top+2*A_bottom);
Izz=2*A_top.*(h-Y_bar).^2+2*A_bottom.*(-Y_bar).^2;
[P,Mz,Mx]=loads(b,c,Sweep_40,D_mx,L,Izz,Z,E);

sig_1=Mz.*(0.95*h-Y_bar)./Izz;
sig_2=Mz.*(h-Y_bar)./Izz;
sig_3=sig_1;
sig_4=Mz.*(0.05*h-Y_bar)./Izz;
sig_5=Mz.*(-Y_bar)./Izz;
sig_6=sig_4;
q=Mx./(2*l.*h);


%-----THIS SECTION COMPUTES THE TOTAL WEIGHT OF A/C-----%

k=6.09375;


%----Point 1----%
T1=P.*(l-a)./l;
tau1_T=T1./(h.*t2);
tau1=q./t2+tau1_T;
sig_eq1=sqrt(sig_1.^2+3*tau1.^2);
sig_cr1=((pi^2)*E*4/(12*(1-nu^2)))*(teq2./(0.95*h)).^2;
tau_cr1=((pi^2)*E*5.5/(12*(1-nu^2)))*(teq2./(0.95*h)).^2;
G(1:3)=k*sig_eq1;
G(4:6)=k*(((sig_1)./sig_cr1)+(tau1./tau_cr1).^2);
G(7:9)=k*((-(sig_1)./sig_cr1)+(tau1./tau_cr1).^2);

%----Point 2----%

tau2=q./t1;
sig_eq2=sqrt(sig_2.^2+3*tau2.^2);
sig_cr2=((pi^2)*E*4/(12*(1-nu^2)))*(teq1./l).^2;
tau_cr2=((pi^2)*E*5.5/(12*(1-nu^2)))*(teq1./l).^2;
G(10:12)=k*sig_eq2;
G(13:15)=k*(((sig_2)./sig_cr2)+(tau2./tau_cr2).^2);
G(16:18)=k*((-(sig_2)./sig_cr2)+(tau2./tau_cr2).^2);

%----Point 3----%

T2=P.*a./l;
tau3_T=-T2./(h.*t2);
tau3=q./t2+tau3_T;
sig_eq3=sqrt(sig_3.^2+3*tau3.^2);
sig_cr3=sig_cr1;
tau_cr3=tau_cr1;
G(19:21)=k*sig_eq3;
G(22:24)=k*(((sig_3)./sig_cr3)+(tau3./tau_cr3).^2);
G(25:27)=k*(((-sig_3)./sig_cr3)+(tau3./tau_cr3).^2);

%----Point 4----%

tau4=-q./t2+tau1_T;
sig_eq4=sqrt(sig_4.^2+3*tau4.^2);
G(28:30)=k*sig_eq4;

%----Point 5----%

tau5=q./t3;
sig_eq5=sqrt(sig_5.^2+3*tau5.^2);
sig_cr5=((pi^2)*E*4/(12*(1-nu^2)))*(teq3./l).^2;
tau_cr5=((pi^2)*E*5.5/(12*(1-nu^2)))*(teq3./l).^2;
G(31:33)=k*sig_eq5;
G(34:36)=k*(((sig_5)./sig_cr5)+(tau5./tau_cr5).^2);
G(37:39)=k*(((-sig_5)./sig_cr5)+(tau5./tau_cr5).^2);

%----Point 6----%

tau6=-q./t2+tau3_T;
sig_eq6=sqrt(sig_6.^2+3*tau6.^2);
G(40:42)=k*sig_eq6;


%-----Constraints-----%

Sig_C=65000*144;
Sig_T=65000*144;

G1(1:3)=(G(1:3))./Sig_C-1;
G1(55:57)=-(G(1:3))./Sig_C-1;
G1(4:9)=G(4:9)-1;

G1(10:12)=(G(10:12))./Sig_C-1;
G1(58:60)=-(G(10:12))./Sig_C-1;
G1(13:18)=G(13:18)-1;

G1(19:21)=(G(19:21))./Sig_C-1;
G1(61:63)=-(G(19:21))./Sig_C-1;
G1(22:27)=G(22:27)-1;

G1(28:30)=(G(28:30))./Sig_T-1;
G1(64:66)=-(G(28:30))./Sig_T-1;

G1(31:33)=(G(31:33))./Sig_T-1;
G1(67:69)=-(G(31:33))./Sig_T-1;
G1(34:39)=G(34:39)-1;

G1(40:42)=(G(40:42))./Sig_T-1; 
G1(70:72)=-(G(40:42))./Sig_T-1; 

G1(43:45)=(1/2)*(ts1+ts3)./h-1;
G1(46:48)=t1./(ts1-.1*t1)-1;
G1(49:51)=t2./(ts2-.1*t2)-1;
G1(52:54)=t3./(ts3-.1*t3)-1;


