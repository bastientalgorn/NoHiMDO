function [Ws , Wf, theta] = SBJ_obj_weight(Z,t,ts,lambda,L)

%----Inputs----%
t=t(:)/12;%convert to feet
ts=ts(:)/12;%convert to feet
C=[500.0 16000.0  4.0  4360.0  0.01375  1.0];

t1=t(1:3); 
t2=t(4:6);  
t3=t(7:9); 
ts1=ts(1:3); 
ts2=ts(4:6); 
ts3=ts(7:9); 
G=4000000*144;
E=10600000*144;
rho_alum=0.1*144;
rho_core=0.1*144/10;
rho_fuel=6.5*7.4805;
Fw_at_t=5;

%----Inputs----%
beta=.9;
[c,c_box,Sweep_40,D_mx,b]=Wing_Mod(Z,lambda);
l=0.6*c_box;
h=Z(1)*(beta*c(1:3))-(0.5)*(ts1+ts3);
A_top=(0.5)*(t1.*l)+(1/6)*(t2.*h);
A_bottom=(0.5)*(t3.*l)+(1/6)*(t2.*h);
Y_bar=h.*(2*A_top)./(2*A_top+2*A_bottom);
Izz=2*A_top.*(h-Y_bar).^2+2*A_bottom.*(-Y_bar).^2;
[P,Mz,Mx,bend_twist,Spanel]=loads(b,c,Sweep_40,D_mx,L,Izz,Z,E);


Phi=(Mx./(4*G*(l.*h).^2)).*(l./t1+2*h./t2+l./t3);
aa=size(bend_twist,2);
twist(1:aa/3)=bend_twist(1:aa/3)+Phi(1)*180/pi;
twist((aa/3)+1:aa*2/3)=bend_twist((aa/3)+1:aa*2/3)+Phi(2)*180/pi;
twist((aa*2/3)+1:aa)=bend_twist((aa*2/3)+1:aa)+Phi(3)*180/pi;
deltaL_divby_q=sum(twist.*Spanel*0.1*2);

%-----THIS SECTION COMPUTES THE TOTAL WEIGHT OF A/C-----%

Wtop_alum=(b/4)*(c(1)+c(4))*mean(t1)*rho_alum;
Wbottom_alum=(b/4)*(c(1)+c(4))*mean(t3)*rho_alum;
Wside_alum=(b/2)*mean(h)*mean(t2)*rho_alum;
Wtop_core=(b/4)*(c(1)+c(4))*mean(ts1-t1)*rho_core;
Wbottom_core=(b/4)*(c(1)+c(4))*mean(ts3-t3)*rho_core;
Wside_core=(b/2)*mean(h)*mean(ts2-t2)*rho_core;
W_wingstruct=Wtop_alum+Wbottom_alum+Wside_alum+Wtop_core+Wbottom_core+Wside_core;
W_fuel_wing=mean(h.*l)*(b/3)*(2)*rho_fuel;
Bh=sqrt(Z(8)*Z(7));
W_ht=3.316*((1+(Fw_at_t/Bh))^-2.0)*((L*C(3)/1000)^0.260)*(Z(7)^0.806);
Wf = C(1) + W_fuel_wing;
Ws = C(2) + W_ht + 2*W_wingstruct;
theta = deltaL_divby_q;




