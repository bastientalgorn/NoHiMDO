function [P,Mz,Mx,bend_twist,Spanel]=loads(b,c,Sweep_40,D_mx,L,Izz,Z,E)

np=9;   %----number of panels per halfspan
n=90;
rn = n/np;

h=(b/2)/n;
x=(0:(b/2-h)/(n-1):b/2-h);
x1=(h:(b/2-h)/(n-1):b/2);

%----Calculate Mx, Mz, and P----%

l=(0:(b/2)/np:(b/2)-(b/2)/np);


c1mc4 = c(1)-c(4);
f_all  =(3*b/10)*sqrt(1-( x.^2)/((b/2)^2));
f1_all =(3*b/10)*sqrt(1-(x1.^2)/((b/2)^2));
C= c(4) + 2*( (b/2- x)/b )*c1mc4;
C1=c(4) + 2*( (b/2-x1)/b )*c1mc4;
A_Tot=(h/4)*(C+C1).*(f_all+f1_all);
Area = sum(reshape(A_Tot,rn,np));
Spanel = (h*(rn)/2)*(C(1:10:n-9)+C(10:10:n));






% cos, tan, and cos^-1 of Sweep
cosSweep = cos(Sweep_40*pi/180);
cosInvSweep = 1/cosSweep;
tanCos2Sweep = tan(Sweep_40*pi/180)*cosSweep*cosSweep;



p=L*Area./sum(Area);
% 
% T = zeros(1,np);
% Mb = zeros(1,np);
% for i=1:np
%   T(i)=sum(p(i:np));
%   Mb(i)=sum( p(i+1:np).*(l(i+1:np)-l(i)) )*cosInvSweep;
% end


% Replace T by:
Tcsp = cumsum(p);
Tsp = Tcsp(end);
T = Tsp-[0 Tcsp(1:end-1)];
pl = p.*l;
Tcspl = cumsum(pl);
Tspl = Tcspl(end);
Mb = ( (Tspl - Tcspl)-l.*(Tsp-Tcsp) )*cosInvSweep;





P=T(1:np/3:np);
P=P(:);
Mx=P.*D_mx;
Mz=Mb(1:np/3:np);
Mz=Mz(:);

%----Calculate Wing Twist due to Bending----%

chord=c(4)+(2*(b/2-l)./b)*c1mc4;
y(1,:)=(l-.4*chord.*tanCos2Sweep)*cosInvSweep;
y(2,:)=(l+.6*chord.*tanCos2Sweep)*cosInvSweep;
y(2,1)=0;
I(1:np/3)=sqrt((Izz(1)^2+Izz(2)^2)/2);
I(np/3+1:2*np/3)=sqrt((Izz(2)^2+Izz(3)^2)/2);
I(2*np/3+1:np)=sqrt((Izz(3)^2)/2);

La=y(1,2:np)-y(1,1:np-1);
La=[0 La];
Lb=y(2,2:np)-y(2,1:np-1);
Lb=[0 Lb];
A=T.*La.^3./(3*E*I)+Mb.*La.^2./(2*E*I);
B=T.*Lb.^3./(3*E*I)+Mb.*Lb.^2./(2*E*I);
Slope_A=T.*La.^2./(2*E*I)+Mb.*La./(E*I);
Slope_B=T.*Lb.^2./(2*E*I)+Mb.*Lb./(E*I);
for i=1:np-1
  Slope_A(i+1)=Slope_A(i)+Slope_A(i+1);
  Slope_B(i+1)=Slope_B(i)+Slope_B(i+1);
  A(i+1)=A(i)+Slope_A(i)*La(i+1)+A(i+1);
  B(i+1)=B(i)+Slope_B(i)*Lb(i+1)+B(i+1);
end
bend_twist=((B-A)./chord)*180/pi;
for i=2:size(bend_twist,2)
   if bend_twist(i)<bend_twist(i-1)
      bend_twist(i)=bend_twist(i-1);
   end
end
%delta_L=bend_twist.*Spanel*q*0.1;
%p_twist=p+delta_L;
%r=1:45/30:45;
%plot(A)
%hold on
%plot(B)
