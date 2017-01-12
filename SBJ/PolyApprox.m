function [FF, Ao, Ai, Aij] = PolyApprox(S, S_new, flag, S_bound)
     for i = 1:size(S,2)
          S_norm(i) = S_new(i)/S(i);  %normalize new S with initial S
          if S_norm(i)>1.25 
               S_norm(i)=1.25;
          elseif S_norm(i)<0.75
               S_norm(i)=0.75;
          end       
          S_shifted(i) = S_norm(i) - 1;  %shift S vector near origin
          
          %DETERMINE BOUNDS ON FF DEPENDING ON SLOPE-SHAPE 
          a=0.1;  
          b=a;
          if flag(i)==5
 
               %CALCULATE POLYNOMIAL COEFFICIENTS (S-ABOUT ORIGIN)
               So=0;
               Sl=So-S_bound(i);
               Su=So+S_bound(i);
               Mtx_shifted = [1 Sl Sl^2; 1 So So^2; 1 Su Su^2];
               F_bound = [1+(.5*a)^2; 1; 1+(.5*b)^2];
               A = Mtx_shifted\F_bound;
               Ao = A(1);
               Ai(i) = A(2);
               Aij(i,i) = A(3);
               %CALCULATE POLYNOMIAL COEFFICIENTS

          else
               switch (flag(i))
                    case 0
                         S_shifted(i) = 0;
                    case 3
                         a=-a;
                         b=a;
                    case 2
                         b=2*a;
                    case 4
                         a=-a;
                         b=2*a;    
               end
               %DETERMINE BOUNDS ON FF DEPENDING ON SLOPE-SHAPE
   
               %CALCULATE POLYNOMIAL COEFFICIENTS (S-ABOUT ORIGIN)
               So=0;
               Sl=So-S_bound(i);
               Su=So+S_bound(i);
               Mtx_shifted = [1 Sl Sl^2; 1 So So^2; 1 Su Su^2];
               F_bound = [1-.5*a; 1; 1+.5*b];
               A = Mtx_shifted\F_bound;
               Ao = A(1);
               Ai(i) = A(2);
               Aij(i,i) = A(3);
               %CALCULATE POLYNOMIAL COEFFICIENTS
          end
     end
 
     %FILL UP OFF DIAGONALS OF Aij 

     R = [...
    0.2736    0.3970    0.8152    0.9230    0.1108
    0.4252    0.4415    0.6357    0.7435    0.1138
    0.0329    0.8856    0.8390    0.3657    0.0019
    0.0878    0.7248    0.1978    0.0200    0.0169
    0.8955    0.4568    0.8075    0.9239    0.2525];

     for i = 1:size(S,2)
          for j = (i+1):size(S,2)
               Aij(i,j) = Aij(i,i)*R(i,j);
               Aij(j,i) = Aij(i,j);
          end
     end 
     
     %CALCULATE POLYNOMIAL
     FF = Ao + Ai*(S_shifted') + (1/2)*(S_shifted)*(Aij)*(S_shifted');    
