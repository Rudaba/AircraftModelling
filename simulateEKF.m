function [pStatesOUT,pCovarianceOUT,innovation,SCovOUT]  = simulateEKF(u,Vt,ALP,BET,accelerations,omega,Height)
%#eml

persistent pStates;
persistent pCovariance;
persistent covDiag;
persistent prevOmegaDots;

if isempty(pStates)
    % set initial state
    pStates = zeros(30,1);
    covDiag = zeros(30,1);
    pStates(1)  = accelerations(1); covDiag(1) = 0.01^2;
    pStates(2)  = accelerations(2); covDiag(2) = 0.01^2;
    pStates(3)  = accelerations(3); covDiag(3) = 0.01^2;
    
    pStates(4)  = omega(1); covDiag(4) = 0.01^2;
    pStates(5)  = omega(2); covDiag(5) = 0.01^2;
    pStates(6)  = omega(3); covDiag(6) = 0.01^2;
    
    
    pStates(7)  = 9.5e-4;   covDiag(7) = 1e-5^2; % CX_dE1 = 9.5e-4;
    pStates(8)  = 8.5e-7;   covDiag(8) = 1e-8^2; % CX_dE2 = 8.5e-7;
    pStates(9)  = 1.75e-4;  covDiag(9) = 1e-5^2; % CY_dE1 = 1.75e-4;
    pStates(10) = 1.55e-3;  covDiag(10) = 1e-4^2; % CY_dR1 = 1.55e-3;
    pStates(11) = 8e-6;     covDiag(11) = 1e-7^2; % CY_dR2 = 8e-6;
    pStates(12) = 4.76e-3;  covDiag(12) = 1e-4^2; % CZ_dE1 = 4.76e-3;
    pStates(13) = 3.3e-5;   covDiag(13)= 1e-6^2; % CZ_dE2 = 3.3e-5;
    pStates(14) = 7.5e-5;   covDiag(14)= 1e-6^2; % CZ_dA1 = 7.5e-5;
    pStates(15) = 6.1e-4;   covDiag(15)= 1e-5^2; % Cl_dA1 = 6.1e-4;
    pStates(16) = 2.5e-5;   covDiag(16)= 1e-6^2; % Cl_dA2 = 2.5e-5;
    pStates(17) = 2.6e-6;   covDiag(17)= 1e-7^2; % Cl_dA3 = 2.6e-6;
    pStates(18) = -2.3e-4;  covDiag(18)= 1e-5^2; % Cl_dR1 = -2.3e-4;
    pStates(19) = 4.5e-6;   covDiag(19)= 1e-7^2; % Cl_dR2 = 4.5e-6;
    pStates(20) = 5.24e-5;  covDiag(20)= 1e-6^2; % Cl_dE1 = 5.24e-5;
    pStates(21) = 6.54e-3;  covDiag(21)= 1e-4^2; % Cm_dE1 = 6.54e-3;
    pStates(22) = 8.49e-5;  covDiag(22)= 1e-6^2; % Cm_dE2 = 8.49e-5;
    pStates(23) = 3.74e-6;  covDiag(23)= 1e-7^2; % Cm_dE3 = 3.74e-6;
    pStates(24) =  3.5e-5;  covDiag(24)= 1e-6^2; % Cm_dA1 = 3.5e-5;
    pStates(25) = 1.4e-5;   covDiag(25)= 1e-6^2; % Cn_dA1 = 1.4e-5;
    pStates(26) = 7.0e-6;   covDiag(26)= 1e-7^2; % Cn_dA2 = 7.0e-6;
    pStates(27) = 8.73e-5;  covDiag(27)= 1e-6^2; % Cn_dE1 = 8.73e-5;
    pStates(28) = 8.7e-6;   covDiag(28)= 1e-7^2; % Cn_dE2 = 8.7e-6;
    pStates(29) = 9.0e-4;   covDiag(29)= 1e-5^2; % Cn_dR1 = 9.0e-4;
    pStates(30) = 4.0e-6;   covDiag(30)= 1e-7^2; % Cn_dR2 = 4.0e-6;
    
    pStates(7:30,1)   = ones(24,1);
    covDiag(7:end,1)  = 1e-2^2;
    %     covDiag(12)       = 1e-5^2;
    
    pCovariance       = diag(covDiag*1e-1);
    
    prevOmegaDots     = [0;0;0];
end


%*****Define Q and R matrices*****
% Qcov = diag(1e-9*[0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2]);
Qcov        = (eye(30,30)*0.1).^2;
Qcov(1,1)   = 0.1^2;
Qcov(2,2)   = 0.1^2;
Qcov(3,3)   = 0.1^2;
Qcov(4,4)   = 0.1^2;
Qcov(5,5)   = 0.1^2;
Qcov(6,6)   = 0.1^2;

Rcov        = diag([0.1^2,0.1^2,0.1^2,0.1^2,0.1^2,0.1^2]);

%*****Constants*****
R2D     = 180/pi;
g       = 9.81;
vs      = 340.3;
R2D     = 180/pi;
rho     = 1.225;%1.225*(1-H/44.331/1000)^4.256;
dt      = 0.005;
m       = 38924 * 0.453592 / 15;

rho     = 1.225;
Ix      = 24970*14.593903*0.092903/15;
Iy      = 122190*14.593903*0.092903/15;
Iz      = 139800*14.593903*0.092903/15;
Ixz     = 1175*14.593903*0.092903/15;
b       = 20/3;
cbar    = 3;
S       = 20;

xCG_ref = 0;
xCG     = 0;

dA          = u(1,1)*180/pi;
dE          = u(2,1)*180/pi;
dR          = u(3,1)*180/pi;

throttle    = u(4,1);

%***Calculate Dynamic Pressure***
qbar = 0.5*rho*Vt^2;

%***Calculate dimensional values***
hT = Height/3048;

Tmax = ((30.21-0.668*hT-6.877*hT^2+1.951*hT^3-0.1512*hT^4) + ...
    (Vt/vs).*(-33.8+3.347*hT+18.13*hT^2-5.865*hT^3+0.4757*hT^4) + ...
    (Vt/vs)^2.*(100.8-77.56*hT+5.441*hT^2+2.864*hT^3-0.3355*hT^4) + ...
    (Vt/vs)^3.*(-78.99+101.4*hT-30.28*hT^2+3.236*hT^3-0.1089*hT^4) + ...
    (Vt/vs)^4.*(18.74-31.6*hT+12.04*hT^2-1.785*hT^3+0.09417*hT^4))*4448.22/20; % Newton's

T = Tmax*throttle;


%*****Extract States From Filter*****
p           = pStates(4);
q           = pStates(5);
r           = pStates(6);

prev_pdot   = prevOmegaDots(1,1);
prev_qdot   = prevOmegaDots(2,1);
prev_rdot   = prevOmegaDots(3,1);

CX_dE1_mp   = pStates(7);
CX_dE2_mp   = pStates(8);
CY_dE1_mp   = pStates(9);
CY_dR1_mp   = pStates(10);
CY_dR2_mp   = pStates(11);
CZ_dE1_mp   = pStates(12);
CZ_dE2_mp   = pStates(13);
CZ_dA1_mp   = pStates(14);
Cl_dA1_mp   = pStates(15);
Cl_dA2_mp   = pStates(16);
Cl_dA3_mp   = pStates(17);
Cl_dR1_mp   = pStates(18);
Cl_dR2_mp   = pStates(19);
Cl_dE1_mp   = pStates(20);
Cm_dE1_mp   = pStates(21);
Cm_dE2_mp   = pStates(22);
Cm_dE3_mp   = pStates(23);
Cm_dA1_mp   = pStates(24);
Cn_dA1_mp   = pStates(25);
Cn_dA2_mp   = pStates(26);
Cn_dE1_mp   = pStates(27);
Cn_dE2_mp   = pStates(28);
Cn_dR1_mp   = pStates(29);
Cn_dR2_mp   = pStates(30);

CX_dE1      = CX_dE1_mp * 9.5e-4;  % CX_dE1 = 9.5e-4;
CX_dE2      = CX_dE2_mp * 8.5e-7;  % CX_dE2 = 8.5e-7;
CY_dE1      = CY_dE1_mp * 1.75e-4; % CY_dE1 = 1.75e-4;
CY_dR1      = CY_dR1_mp * 1.55e-3; % CY_dR1 = 1.55e-3;
CY_dR2      = CY_dR2_mp * 8e-6;    % CY_dR2 = 8e-6;
CZ_dE1      = CZ_dE1_mp * 4.76e-3; % CZ_dE1 = 4.76e-3;
CZ_dE2      = CZ_dE2_mp * 3.3e-5;  % CZ_dE2 = 3.3e-5;
CZ_dA1      = CZ_dA1_mp * 7.5e-5;  % CZ_dA1 = 7.5e-5;
Cl_dA1      = Cl_dA1_mp * 6.1e-4;  % Cl_dA1 = 6.1e-4;
Cl_dA2      = Cl_dA2_mp * 2.5e-5;  % Cl_dA2 = 2.5e-5;
Cl_dA3      = Cl_dA3_mp * 2.6e-6;  % Cl_dA3 = 2.6e-6;
Cl_dR1      = Cl_dR1_mp * -2.3e-4; % Cl_dR1 = -2.3e-4;
Cl_dR2      = Cl_dR2_mp * 4.5e-6;  % Cl_dR2 = 4.5e-6;
Cl_dE1      = Cl_dE1_mp * 5.24e-5; % Cl_dE1 = 5.24e-5;
Cm_dE1      = Cm_dE1_mp * 6.54e-3; % Cm_dE1 = 6.54e-3;
Cm_dE2      = Cm_dE2_mp * 8.49e-5; % Cm_dE2 = 8.49e-5;
Cm_dE3      = Cm_dE3_mp * 3.74e-6; % Cm_dE3 = 3.74e-6;
Cm_dA1      = Cm_dA1_mp * 3.5e-5;  % Cm_dA1 = 3.5e-5;
Cn_dA1      = Cn_dA1_mp * 1.4e-5;  % Cn_dA1 = 1.4e-5;
Cn_dA2      = Cn_dA2_mp * 7.0e-6;  % Cn_dA2 = 7.0e-6;
Cn_dE1      = Cn_dE1_mp * 8.73e-5; % Cn_dE1 = 8.73e-5;
Cn_dE2      = Cn_dE2_mp * 8.7e-6;  % Cn_dE2 = 8.7e-6;
Cn_dR1      = Cn_dR1_mp * 9.0e-4;  % Cn_dR1 = 9.0e-4;
Cn_dR2      = Cn_dR2_mp * 4.0e-6;  % Cn_dR2 = 4.0e-6;


% CX = -0.0434 + 2.39e-3*ALP+2.53e-5*BET^2-1.07e-6*ALP*BET^2+9.5e-4*dE-8.5e-7*dE*BET^2+(180*q*cbar/pi/2/Vt)*(8.73e-3+0.001*ALP-1.75e-4*ALP^2);
% CY = -0.012*BET+1.55e-3*dR-8e-6*dR*ALP+(180*b/pi/2/Vt)*(2.25e-3*p+0.0117*r-3.67e-4*r*ALP+1.75e-4*r.*dE);
% CZ = -0.131-0.0538*ALP-4.76e-3*dE-3.3e-5*dE*ALP-7.5e-5*dA.^2+(180*q*cbar/pi/2/Vt).*(-0.111+5.17e-3*ALP-1.1e-3*ALP^2);
% Cl = -5.98e-4*BET-2.83e-4*ALP*BET+1.51e-5*ALP^2*BET-dA.*(6.1e-4+2.5e-5*ALP-2.6e-6*ALP^2)-dR.*(-2.3e-4+4.5e-6*ALP)+(180*b/2/pi/Vt).*(-4.12e-3*p-5.24e-4*p*ALP+4.36e-5*p*ALP^2+4.36e-4*r+1.05e-4*r*ALP+5.24e-5*r.*dE);
% Cm = -6.61e-3-2.67e-3*ALP-6.48e-5*BET^2-2.65e-6*ALP*BET^2-6.54e-3*dE-8.49e-5*dE*ALP+3.74e-6*dE*BET^2-3.5e-5*dA.^2+(180*q*cbar/2/pi/Vt).*(-0.0473-1.57e-3*ALP)+(xCG_ref-xCG)*CZ;
% Cn = 2.28e-3*BET+1.79e-6*BET^3+1.4e-5*dA+7.0e-6*dA*ALP-9.0e-4*dR+4.0e-6*dR*ALP+(180*b/2/pi/Vt).*(-6.63e-5*p-1.92e-5*p*ALP+5.06e-6*p*ALP^2-6.06e-3*r-8.73e-5*r.*dE+8.7e-6*r.*dE*ALP)-cbar/b*(xCG_ref-xCG)*CY;

CX = -0.0434 + 2.39e-3*ALP+2.53e-5*BET^2-1.07e-6*ALP*BET^2+CX_dE1*dE-CX_dE2*dE*BET^2+(180*q*cbar/pi/2/Vt)*(8.73e-3+0.001*ALP-1.75e-4*ALP^2);
CY = -0.012*BET+CY_dR1*dR-CY_dR2*dR*ALP+(180*b/pi/2/Vt)*(2.25e-3*p+0.0117*r-3.67e-4*r*ALP+CY_dE1*r.*dE);
CZ = -0.131-0.0538*ALP-CZ_dE1*dE-CZ_dE2*dE*ALP-CZ_dA1*dA.^2+(180*q*cbar/pi/2/Vt).*(-0.111+5.17e-3*ALP-1.1e-3*ALP^2);
Cl = -5.98e-4*BET-2.83e-4*ALP*BET+1.51e-5*ALP^2*BET-dA.*(Cl_dA1+Cl_dA2*ALP-Cl_dA3*ALP^2)-dR.*(Cl_dR1+Cl_dR2*ALP)+(180*b/2/pi/Vt).*(-4.12e-3*p-5.24e-4*p*ALP+4.36e-5*p*ALP^2+4.36e-4*r+1.05e-4*r*ALP+Cl_dE1*r.*dE);
Cm = -6.61e-3-2.67e-3*ALP-6.48e-5*BET^2-2.65e-6*ALP*BET^2-Cm_dE1*dE-Cm_dE2*dE*ALP+Cm_dE3*dE*BET^2-Cm_dA1*dA.^2+(180*q*cbar/2/pi/Vt).*(-0.0473-1.57e-3*ALP)+(xCG_ref-xCG)*CZ;
Cn = 2.28e-3*BET+1.79e-6*BET^3+Cn_dA1*dA+Cn_dA2*dA*ALP-Cn_dR1*dR+Cn_dR2*dR*ALP+(180*b/2/pi/Vt).*(-6.63e-5*p-1.92e-5*p*ALP+5.06e-6*p*ALP^2-6.06e-3*r-Cn_dE1*r.*dE+Cn_dE2*r.*dE*ALP)-cbar/b*(xCG_ref-xCG)*CY;

c1 = ((Iy-Iz)*Iz-Ixz^2)/(Ix*Iz-Ixz^2);
c2 = (Ix-Iy+Iz)*Ixz/(Ix*Iz-Ixz^2);
c3 = Iz/(Ix*Iz-Ixz^2);
c4 = Ixz/(Ix*Iz-Ixz^2);
c5 = (Iz-Ix)/Iy;
c6 = Ixz/Iy;
c7 = 1/Iy;
c8 = ((Ix-Iy)*Ix-Ixz^2)/(Ix*Iz-Ixz^2);
c9 = Ix/(Ix*Iz-Ixz^2);

heng = 0;
pdot = (c1*r+c2*p+c4*heng).*q+qbar*S*b*(c3*Cl+c4*Cn);
qdot = (c5*p-c7*heng).*r-c6*(p.^2-r.^2)+qbar*S*cbar*c7*Cm;
rdot = (c8*p-c2*r+c9*heng).*q+qbar*S*b*(c4*Cl+c9*Cn);

%**************Predict************
%Accelerations
pStates(1,1)    = (qbar*S*CX+T)/m; 
pStates(2,1)    = qbar*S*CY/m; 
pStates(3,1)    = qbar*S*CZ/m;

%Angular Rates
%Trapezoidal Integration and Jacobian
pStates(4,1)    = p + dt*(prev_pdot+pdot)/2;
pStates(5,1)    = q + dt*(prev_qdot+qdot)/2;
pStates(6,1)    = r + dt*(prev_rdot+rdot)/2;
prevOmegaDots   = [pdot;qdot;rdot];

F = eye(30);

F(1,5)  = (1.460334333819828426948389041052*qbar*(- 0.000175*ALP^2 + 0.001*ALP + 0.00873))/Vt;
F(1,7)  = 0.016991761536715991584456999092618*dE*qbar;
F(1,8)  = -0.016991761536715991584456999092618*BET^2*dE*qbar;

F(2,4)  = (0.0073016716690991420258687569821225*qbar)/Vt;
F(2,6)  = (3.2451874084885075670527808809433*qbar*(CY_dE1*dE - 0.00036699999999999997598101875162513*ALP + 0.0117))/Vt;
F(2,9)  = 0.016991761536715991584456999092618*dR*qbar;
F(2,10) = -0.016991761536715991584456999092618*ALP*dR*qbar;
F(2,11) = (3.2451874084885075670527808809433*dE*qbar*r)/Vt;

F(3,5)  = -(1.460334333819828426948389041052*qbar*(0.0011*ALP^2 - 0.00517*ALP + 0.111))/Vt;
F(3,12) = -0.016991761536715991584456999092618*dE*qbar;
F(3,13) = -0.016991761536715991584456999092618*ALP*dE*qbar;
F(3,14) = -0.016991761536715991584456999092618*dA^2*qbar;

F(4,4)  = 0.5*dt*(0.014338034095370216433606991301986*q - 133.33333333333333333333333333333*qbar*((0.084653457277151632633481163237447*(- 0.000043600000000000002685698885507293*ALP^2 + 0.00052400000000000005167394290239713*ALP + 0.00412))/Vt + (0.00071150080329508703831251117964528*(- 0.000005059999999999999834552916883057*ALP^2 + 0.000019199999999999999131163747057016*ALP + 0.000066299999999999998799744826971647))/Vt)) + 1.0;
F(4,5)  = 0.5*dt*(0.014338034095370216433606991301986*p - 0.705920992793835466727614402771*r);
F(4,6)  = -0.5*dt*(0.705920992793835466727614402771*q + 133.33333333333333333333333333333*qbar*((0.00071150080329508703831251117964528*(Cn_dE1*dE - 1.0*ALP*Cn_dE2*dE + 0.00606))/Vt - (0.084653457277151632633481163237447*(0.00010500000000000000435415592470179*ALP + Cl_dE1*dE + 0.00043600000000000002685698885507293))/Vt));
F(4,15) = -0.029549631053652999392034050885059*dA*dt*qbar;
F(4,16) = -0.029549631053652999392034050885059*ALP*dA*dt*qbar;
F(4,17) = 0.029549631053652999392034050885059*ALP^2*dA*dt*qbar;
F(4,18) = -0.029549631053652999392034050885059*dR*dt*qbar;
F(4,19) = -0.029549631053652999392034050885059*ALP*dR*dt*qbar;
F(4,20) = (5.6435638184767755088987442158298*dE*dt*qbar*r)/Vt;
F(4,25) = 0.00024836063296167578075904676844961*dA*dt*qbar;
F(4,26) = 0.00024836063296167578075904676844961*ALP*dA*dt*qbar;
F(4,27) = -0.00024836063296167578075904676844961*dR*dt*qbar;
F(4,28) = 0.00024836063296167578075904676844961*ALP*dR*dt*qbar;
F(4,29) = -(0.047433386886339135887500745309686*dE*dt*qbar*r)/Vt;
F(4,30) = (0.047433386886339135887500745309686*ALP*dE*dt*qbar*r)/Vt;

F(5,4)  = -0.5*dt*(0.01923234307226450609706195269662*p - 0.93976593829282265324494639495867*r);
F(5,5)  = 1.0 - (0.23344767165356548584164144954176*dt*qbar*(0.00157*ALP + 0.0473))/Vt;
F(5,6)  = 0.5*dt*(0.93976593829282265324494639495867*p + 0.01923234307226450609706195269662*r);
F(5,21) = -0.0027162870009795687031503574893065*dE*dt*qbar;
F(5,22) = -0.0027162870009795687031503574893065*ALP*dE*dt*qbar;
F(5,23) = 0.0027162870009795687031503574893065*BET^2*dE*dt*qbar;
F(5,24) = -0.0027162870009795687031503574893065*dA^2*dt*qbar;

F(6,4)  = -0.5*dt*(0.69609284164731566324491041086731*q + 133.33333333333333333333333333333*qbar*((0.00071150080329508703831251117964528*(- 0.000043600000000000002685698885507293*ALP^2 + 0.00052400000000000005167394290239713*ALP + 0.00412))/Vt + (0.015120148985768787881099866668956*(- 0.000005059999999999999834552916883057*ALP^2 + 0.000019199999999999999131163747057016*ALP + 0.000066299999999999998799744826971647))/Vt));
F(6,5)  = -0.5*dt*(0.69609284164731566324491041086731*p + 0.014338034095370216433606991301986*r);
F(6,6)  = 1.0 - 0.5*dt*(0.014338034095370216433606991301986*q + 133.33333333333333333333333333333*qbar*((0.015120148985768787881099866668956*(Cn_dE1*dE - 1.0*ALP*Cn_dE2*dE + 0.00606))/Vt - (0.00071150080329508703831251117964528*(0.00010500000000000000435415592470179*ALP + Cl_dE1*dE + 0.00043600000000000002685698885507293))/Vt));
F(6,15) = -0.00024836063296167578075904676844961*dA*dt*qbar;
F(6,16) = -0.00024836063296167578075904676844961*ALP*dA*dt*qbar;
F(6,17) = 0.00024836063296167578075904676844961*ALP^2*dA*dt*qbar;
F(6,18) = -0.00024836063296167578075904676844961*dR*dt*qbar;
F(6,19) = -0.00024836063296167578075904676844961*ALP*dR*dt*qbar;
F(6,20) = (0.047433386886339135887500745309686*dE*dt*qbar*r)/Vt;
F(6,25) = 0.0052779276638749319898458178812461*dA*dt*qbar;
F(6,26) = 0.0052779276638749319898458178812461*ALP*dA*dt*qbar;
F(6,27) = -0.0052779276638749319898458178812461*dR*dt*qbar;
F(6,28) = 0.0052779276638749319898458178812461*ALP*dR*dt*qbar;
F(6,29) = -(1.0080099323845858587399911112638*dE*dt*qbar*r)/Vt;
F(6,30) = (1.0080099323845858587399911112638*ALP*dE*dt*qbar*r)/Vt;


%*****Predict Covariance Forward*****
pCovariance = F*pCovariance*F'+Qcov;

%*****Update*****
H               = eye(6,30);

innovation      = [accelerations;omega] - pStates(1:6,1);
SCov            = H*pCovariance*H'+ Rcov;
k               = pCovariance*H'*inv(SCov);
correction      = k*innovation;
pStates         = pStates + correction;
pCovariance     = pCovariance - k*H*pCovariance;


pCovarianceOUT  = diag(pCovariance);
pStatesOUT      = pStates;
SCovOUT         = diag(SCov); 