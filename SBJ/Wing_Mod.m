function[c,c_box,Sweep_40,D_mx,b,l]=Wing_Mod(Z,lambda)

b=max(2,real(sqrt(Z(4)*Z(6))));
c(1)=2*Z(6)/((1+lambda)*b);
c(4)=lambda*c(1);
x(1)=0;
y(1)=0;
x(2)=c(1);
y(2)=0;
x(7)=(b/2)*tan(Z(5)*pi/180);
y(7)=b/2;
x(8)=x(7)+c(4);
y(8)=b/2;
y(3)=b/6;
x(3)=(x(7)/y(7))*y(3);
y(5)=b/3;
x(5)=(x(7)/y(7))*y(5);
x(6)=x(8)+((x(2)-x(8))/y(8))*(y(8)-y(5));
y(6)=y(5);
x(4)=x(8)+((x(2)-x(8))/y(8))*(y(8)-y(3));
y(4)=y(3);
c(2)=x(4)-x(3);
c(3)=x(6)-x(5);
TE_sweep=(atan((x(8)-x(2))/y(8)))*180/pi;
Sweep_40=(atan(((x(8)-0.6*(x(8)-x(7)))-0.4*x(2))/y(8)))*180/pi;

l=(.4*c(1:3))*cos(Z(5)*pi/180);
k=(.6*c(1:3)*sin((90-TE_sweep)*pi/180)/sin((90+TE_sweep-Z(5))*pi/180));
c_box=l+k;
D_mx=l-0.407*c_box;

l=l(:); k=k(:); D_mx=D_mx(:); c=c(:); c_box=c_box(:);
