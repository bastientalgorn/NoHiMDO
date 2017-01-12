function G3 = SBJ_constraint_power(h,M,T,D)

Dim_Throttle = T*16168;   %--non-diminsional throttle setting

%-----THIS SECTION COMPUTES POLYNOMIAL CONSTRAINT FUNCTIONS-----%

S_initial1=[M,h,D];
S1=[M,h,T];
flag1 = [2,4,2];
bound1 = [.25,.25,.25];
Temp_uA=1.02;
G3(1)=PolyApprox(S_initial1,S1,flag1,bound1) /Temp_uA-1; %--engine temperature


p=[11483.7822254806 10856.2163466548 -0.5080237941 3200.157926969 -0.1466251679 0.0000068572];  
Throttle_uA=p(1)+p(2)*M+p(3)*h+p(4)*M^2+2*p(5)*M*h+p(6)*h^2;
G3(2)=Dim_Throttle/Throttle_uA-1;  %--throttle setting
