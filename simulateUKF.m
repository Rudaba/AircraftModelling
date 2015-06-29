function [pStatesOUT,pCovarianceOUT,innovation,SCovOUT]  = simulateUKF(u,Vt,ALP,BET,accelerations,omega,Height)
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
%Calculate sigma points
kapa            = 0.1;
na              = length(pStates); %length of states
N               = size(R,1); %length of measurements
sigma           = zeros(na,2*na+1);
sigma_meas      = zeros(N,2*na+1);
sigma(:,1)      = pStates;
sigma_meas(:,1) = X_Filter(1:3);
W               = zeros(1+2*na,1);
W(1,1)          = kapa/(na+kapa);

for i = 1:2*length(pStates)
    eqn = chol((na+kapa)*pCovariance);
    if i>=1 && i<= na
        sigma(:,i+1) = pStates + eqn(i,:)';
    else
        sigma(:,i+1) = pStates - eqn(i-na,:)';
    end
    W(i+1) = 1/(2*(na+kapa));
end

%Time update for each sigma point
for i = 1:2*na+1
    
    p           = sigma(4,i);
    q           = sigma(5,i);
    r           = sigma(6,i);
    
    prev_pdot   = prevOmegaDots(1,1);
    prev_qdot   = prevOmegaDots(2,1);
    prev_rdot   = prevOmegaDots(3,1);
    
    CX_dE1_mp   = sigma(7,i);
    CX_dE2_mp   = sigma(8,i);
    CY_dE1_mp   = sigma(9,i);
    CY_dR1_mp   = sigma(10,i);
    CY_dR2_mp   = sigma(11,i);
    CZ_dE1_mp   = sigma(12,i);
    CZ_dE2_mp   = sigma(13,i);
    CZ_dA1_mp   = sigma(14,i);
    Cl_dA1_mp   = sigma(15,i);
    Cl_dA2_mp   = sigma(16,i);
    Cl_dA3_mp   = sigma(17,i);
    Cl_dR1_mp   = sigma(18,i);
    Cl_dR2_mp   = sigma(19,i);
    Cl_dE1_mp   = sigma(20,i);
    Cm_dE1_mp   = sigma(21,i);
    Cm_dE2_mp   = sigma(22,i);
    Cm_dE3_mp   = sigma(23,i);
    Cm_dA1_mp   = sigma(24,i);
    Cn_dA1_mp   = sigma(25,i);
    Cn_dA2_mp   = sigma(26,i);
    Cn_dE1_mp   = sigma(27,i);
    Cn_dE2_mp   = sigma(28,i);
    Cn_dR1_mp   = sigma(29,i);
    Cn_dR2_mp   = sigma(30,i);
    
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
    sigma(1,i)    = (qbar*S*CX+T)/m;
    sigma(2,i)    = qbar*S*CY/m;
    sigma(3,i)    = qbar*S*CZ/m;
    
    %Angular Rates
    %Trapezoidal Integration and Jacobian
    sigma(4,i)    = p + dt*(prev_pdot+pdot)/2;
    sigma(5,i)    = q + dt*(prev_qdot+qdot)/2;
    sigma(6,i)    = r + dt*(prev_rdot+rdot)/2;
    prevOmegaDots   = [pdot;qdot;rdot];
end



%Time update for each sigma point
for i = 1:2*na+1
    
    x   = sigma(1,i);
    y   = sigma(2,i);
    psi = sigma(3,i);
    RR  = sigma(4,i);
    RL  = sigma(5,i);
    
    psiDot          = RR/(2*b)*omegaR - RL/(2*b)*omegaL;
    sigma(3,i)      = psi + psiDot*dt;
    
    V               = RR/2*omegaR + RL/2*omegaL;
    
    xDot            =  V*cos(psi);
    sigma(1,i)      = x + xDot*dt;
    
    
    yDot            = V*sin(psi);
    sigma(2,i)      = y + yDot*dt;
    
    %Measurement predictions
    sigma_meas(1,i) = sigma(1,i);
    sigma_meas(2,i) = sigma(2,i);
    sigma_meas(3,i) = sigma(3,i);
    
end


%Predict states and Convariance forward
pStates     = zeros(na,1);
P           = zeros(na,na);
ymeas       = zeros(N,1);

for i = 1:2*na+1
    pStates = pStates + W(i)*sigma(:,i);
end

for i = 1:2*na+1
    P = P + W(i)*[sigma(:,i)-pStates]*[sigma(:,i)-pStates]';
end

P = P + Q;

for i = 1:2*na+1
    ymeas = ymeas + W(i)*sigma_meas(:,i);
end

%Update
Pxz = zeros(na,N);
Pzz = zeros(N,N);

for i = 1:2*na+1
    Pxz = Pxz + W(i)*[sigma(:,i)-pStates]*[sigma_meas(:,i)-ymeas]';
    Pzz = Pzz + W(i)*[sigma_meas(:,i)-ymeas]*[sigma_meas(:,i)-ymeas]';
end

SCovOUT     = R_Noise + Pzz;
k           = Pxz*inv(S);
innovation  = measurement - ymeas;
correction  = k*innovation;
X_Filter    = X_Filter + correction;
P           = P - k*S*k';