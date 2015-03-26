%Constants
g = 9.81;

%Initial parameters
x=[         0.159912470544201
         0.622926956151104
       -0.0741106103502124
     2.54142987524092e-023
     1.95326498790785e-023];

y0 = [0;0;-100;50;0;0;0;x(1);0;0;0;0]; %N, E, D, Vn, Ve, Vd, phi, theta, psi, P, Q, R
u0 = [x(2);x(3);x(4);x(5)];
% u  = u0; 
% y  = y0; 

%Aircraft Parameters

% omegaDemands    = load('OMEGA_DEMAND');
% 
% omegaTime       = omegaDemands.OD.Data(:,1);
% 
% rollRateDemand  = omegaDemands.OD.Data(:,2);
% 
% pitchRateDemand = omegaDemands.OD.Data(:,3);



