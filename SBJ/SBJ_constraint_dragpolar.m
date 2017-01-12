function G2=SBJ_constraint_dragpolar(Z,Lw,Lh,Wt,tc)

%----Inputs----%
C=[500.0 16000.0  4.0  4360.0  0.01375  1.0];
S_ht=Z(7);     
Nh=C(6);
       
%-----Drag computations----%

if Z(2)<36089
     V = Z(3)*(1116.39*sqrt(1-(6.875e-06*Z(2))));
     rho = (2.377e-03)*(1-(6.875e-06*Z(2)))^4.2561;
else
     V = Z(3)*968.1;
     rho = (2.377e-03)*(.2971)*exp(-(Z(2)-36089)/20806.7);
end
q=.5*rho*(V^2);
%%%% Modified by S. Tosserams:
%%%% scale coefficients for proper conditioning of matrix A 
a=q*Z(6)/1e5;
b=Nh*q*S_ht/1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=Lw;
d=(Lh)*Nh*(S_ht/Z(6));
A=[a b; c d];
%%%% Modified by S. Tosserams:
%%%% scale coefficient Wt for proper conditioning of matrix A 
B=[Wt/1e5; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLo=(A\B);


%-----THIS SECTION COMPUTES CONSTRAINT POLYNOMIALS-----%

S_initial2=tc;
S2=[Z(1)];
flag2=[1];
bound2=[.25];
G(1)=PolyApprox(S_initial2,S2,flag2,bound2);   %--adverse pressure gradient

%------Constraints------%

Pg_uA=1.1;
G2(1)=G(1)/Pg_uA-1;
if CLo(1) > 0
    G2(2)=(2*(CLo(2)))-(CLo(1));
    G2(3)=(2*(-CLo(2)))-(CLo(1));
else
    G2(2)=(2*(-CLo(2)))-(CLo(1));
    G2(3)=(2*(CLo(2)))-(CLo(1));
end





